! ---------------------------------------------------------------------------
! GPU-offloaded 4PCF query using OpenACC.
!
! Two entry points mirror the CPU module:
!   query_graph_4pcf_gpu        — all configurations, no parity channel
!   query_graph_4pcf_parity_gpu — all configurations + parity decomposition
!
! query_graph_4pcf_gpu uses a two-phase "local-matrix" (lmat) strategy:
!
!   Phase 1 — precompute per-hub connectivity matrix (C(nn_i,2) searches):
!     lmat(k2, k1) = distance bin of edge between hub-neighbors k1 < k2.
!     Stored transposed (first index varies) so vector lanes reading
!     lmat(k3, k1) and lmat(k3, k2) with consecutive k3 are coalesced.
!
!   Phase 2 — triple loop using O(1) lmat lookups:
!     Replaces C(nn_i,3)×3 ≈ 162 000 random binary searches (avg nn≈69)
!     with the C(nn_i,2) ≈ 2 346 from Phase 1 → ~69× fewer cache misses.
!
! Parallel mapping: hubs are dealt round-robin onto NGANG_4PCF gangs
! (i = ig, ig+NGANG, ...); within a gang the innermost loop (k2 in Phase 1,
! k3 in Phase 2) runs across VECTOR lanes.  Every lane of a warp then reads
! the same adjacency list / consecutive lmat entries — coherent and coalesced.
! (A bare GANG loop over hubs would let nvfortran put a different hub on each
! vector lane: divergent accesses; with VECTOR_LENGTH(1) the lanes are wasted.)
!
! Memory strategy: lmat must be per-gang.  PRIVATE arrays on a gang loop with
! VECTOR_LENGTH > 1 are allocated PER LANE by nvfortran (→ ILLEGAL_ADDRESS for
! 360 KB × 64 lanes), so instead each gang owns a slice of a global scratch
! array lmat_g(:,:,ig).  The explicit ig loop provides the slot index that
! OpenACC otherwise does not expose.  Accumulation goes to per-gang partials
! part_*(config, ig) via ATOMIC UPDATE (lanes within a gang race; gangs do
! not), summed on the host afterwards.
!
! Runtime sizing: lmat extent = longest CSR row (csr_max_row_len); the gang
! count and edge-window size are derived from free device memory.
!
! Chunked mode (graphs larger than device memory): Phase 2 reads lmat columns
! built from TWO different neighbor rows (id1's and id2's), so tiles run over
! PAIRS of search windows (cw1 <= cw2; id2 > id1 in the sorted rows guarantees
! window(id2) >= window(id1)).  Per tile, Phase 1a builds lmat columns whose
! neighbor id lies in cw1 (searching staged window 1), Phase 1b those in cw2
! (staged window 2), and Phase 2 processes exactly the (id1 in cw1, id2 in cw2)
! pairs — each quadruplet is counted in exactly one tile.  Phase 1 columns are
! rebuilt across tiles (O(nwin) redundancy) but Phase 1 is O(nn^2) against
! Phase 2's O(nn^3), so the overhead is small.  Device footprint: hub window +
! two staged search windows + lmat scratch.
!
! Jackknife (-njk): per-region touching sums N4jk/R4jk (n_configs, nch, njk)
! are split into two parts (see jk_gang_layout):
!
!   * the HUB-region term.  Hubs are permuted so that each region's hubs are
!     contiguous and the gangs are split into one contiguous group per
!     region, so a gang only ever processes hubs of ONE region.  Its
!     per-gang partials part_*(:, ig) are then exactly that region's hub
!     term and are added into N4jk/R4jk on the host after the kernel — no
!     device atomic at all for the (dominant) hub-region contribution.
!
!   * the NEIGHBOUR-region terms: regions of the three neighbours that
!     differ from the hub's region get direct ATOMIC UPDATEs into the
!     device copies of N4jk/R4jk.  Only hubs whose adjacency list crosses
!     a region boundary ("mixed" hubs, flagged in Phase 1) run this code;
!     interior hubs skip the jackknife block entirely.  No slot or gang
!     dimension is used for these: at 4PCF configuration counts that
!     layout would need terabytes, and boundary quadruplets are rare.
!
! The arrays ride the outermost DATA region as COPY, so the accumulated
! sums land straight in the host arrays that write_4pcf_jackknife consumes.
!
! Write routines and S4 symmetry table are reused from query_4pcf_module.
! ---------------------------------------------------------------------------
module query_4pcf_gpu_module
  use kdtree2_precision_module
  use iso_fortran_env, only: int8, int64
  use config_module
  use csr_module
  use query_4pcf_module, only: &
    finish_4pcf_output, &
    dir_x, dir_y, dir_z, VOL_DEGEN_TOL, px, py, pz
  implicit none

contains

  ! ---------------------------------------------------------------------------
  ! 4PCF without parity — gang-slot lmat precomputation + vectorized loops.
  ! ---------------------------------------------------------------------------
  subroutine query_graph_4pcf_gpu(istart, iend)
    integer, intent(in) :: istart, iend

    integer :: ig, i, k1, k2, k3, nn_i, id1, id2, id3, config_idx
    integer(int64) :: base_i
    integer(int8) :: ind1, ind2, ind3, ind4, ind5, ind6
    real(kdkind) :: w12, w4
    logical :: rand12
    integer :: num_data_g, n_configs_g
    ! Jackknife: hub/neighbor region labels (jr1..jr3 hoisted per loop
    ! level and gang-private; jr4 is read in the vector loop).  mixed: some
    ! neighbour of the hub lies in another region; dojk: run the in-kernel
    ! jackknife block for this hub; do2/do3/do4: neighbour region is new.
    integer :: njk_g, jr1, jr2, jr3, jr4
    logical :: mixed, dojk, do2, do3, do4, direct_hub
    ! Gang/region layout (jk_gang_layout): p indexes hub_perm, g the group.
    integer :: p, g
    integer, allocatable :: hub_perm(:), reg_lo(:), gang_grp(:), gang_lo(:), gang_n(:)
    ! Hub-region jackknife terms folded from the per-gang partials on the
    ! host (kept apart from N4jk/R4jk, whose device copies collect the
    ! neighbour-region atomics and are copied back at the end of the DATA
    ! region); folded = partials already folded per hub window (chunked).
    real(kdkind), allocatable :: hub_n4jk(:,:,:), hub_r4jk(:,:,:)
    logical :: folded

    ! Runtime sizing
    integer :: max_nn, ngang, nwin, env_stat
    integer(int64) :: mnn2, lmat_bytes, freeb, win_edges, jk_bytes
    character(32) :: env

    ! Chunking bookkeeping
    integer :: cw1, cw2, hw, do1b
    integer :: w1lo, w1hi, w2lo, w2hi, hublo, hubhi
    integer(int64) :: off1, off2, n1, n2, maxw, elo_h, ehi_h
    integer, allocatable :: split_id(:)
    integer(int64), allocatable :: split_e(:)
    integer, allocatable :: stage1_id(:), stage2_id(:)
    integer(int8), allocatable :: stage1_dist(:), stage2_dist(:)

    ! Per-gang connectivity-matrix scratch; lmat_g(k2, k1, ig) = distance bin
    ! of the edge between hub-neighbors k1 and k2 (k2 > k1).  Transposed so
    ! vector lanes (consecutive k2 / k3) read stride-1 → coalesced.
    integer(int8), allocatable :: lmat_g(:,:,:)
    ! Per-gang partial accumulators, summed over gangs on the host.
    real(kdkind), allocatable :: part_n4(:,:), part_r4(:,:)

    if (.not. cfg%half_graph) then
      print *, 'ERROR: GPU 4PCF requires half_graph=.true.'
      stop 1
    end if

    if (cfg%rank == 0) print *, 'Performing 4PCF (all configs, GPU lmat)'
    if (cfg%rank == 0) print *, 'begin querying the graph'

    num_data_g  = cfg%num_data
    n_configs_g = cfg%n_configs_4pcf
    njk_g       = cfg%njk

    ! ---- Runtime sizing -------------------------------------------------
    max_nn = max(csr_max_row_len(), 1)
    mnn2   = int(max_nn, int64)**2
    freeb  = gpu_free_mem_bytes()

    ! Jackknife accumulators resident on device for the whole query
    ! (N4jk + R4jk, one channel).  Charged against the gang and window
    ! budgets below so legitimate sizes shrink the windows instead of
    ! exhausting the device.
    jk_bytes = 0_int64
    if (njk_g > 0) then
      jk_bytes = 2_int64 * int(n_configs_g, int64) * int(njk_g, int64) * 8_int64
      if (cfg%rank == 0) &
        print '(" 4PCF GPU jackknife accumulators: ",f8.2," GB on device")', &
          real(jk_bytes, kdkind) / 1.0d9
      if (freeb >= 0 .and. jk_bytes > max(freeb - RESERVE_BYTES, 0_int64) / 2) then
        print *, 'ERROR: 4PCF jackknife accumulators exceed half the free'
        print *, '       device memory; reduce -njk or -nbins'
        stop 1
      end if
    end if

    if (freeb < 0) then
      ngang = 1024        ! host fallback: no device limit, modest scratch
    else
      ! lmat scratch capped at ~1/3 of available device memory.
      ngang = int(min(4096_int64, max(512_int64, &
                      ((freeb - RESERVE_BYTES - jk_bytes) / 3) / mnn2)))
    end if
    lmat_bytes = int(ngang, int64) * mnn2

    ! Three edge windows (hub + 2 search) resident in chunked mode.
    win_edges = csr_edge_window(3, lmat_bytes + jk_bytes)
    nwin = int((csr_total_edges - 1) / win_edges) + 1
    ! Without an explicit override, prefer single-pass whenever the whole
    ! CSR fits alongside the scratch (= up to 3 windows' worth of edges).
    call get_environment_variable('GRAMSCI_GPU_WIN_EDGES', env, status=env_stat)
    if (env_stat /= 0 .and. nwin > 1 .and. nwin <= 3) nwin = 1

    if (cfg%rank == 0) &
      print '("4PCF GPU: max_nn=",i0,"  gangs=",i0,"  lmat scratch=",f6.2," GB",'// &
            '"  windows=",i0)', max_nn, ngang, real(lmat_bytes, kdkind)/1.0d9, nwin

    call jk_gang_layout(istart, iend, ngang, njk_g, hub_perm, reg_lo, &
                        gang_grp, gang_lo, gang_n, direct_hub)
    allocate(hub_n4jk(n_configs_g, 1, max(njk_g, 1)))
    allocate(hub_r4jk(n_configs_g, 1, max(njk_g, 1)))
    hub_n4jk = 0.0d0
    hub_r4jk = 0.0d0
    folded = .false.

    allocate(lmat_g(max_nn, max_nn, ngang))
    allocate(part_n4(n_configs_g, ngang))
    allocate(part_r4(n_configs_g, ngang))
    part_n4 = 0.0d0
    part_r4 = 0.0d0

    if (nwin == 1) then
      ! =====================================================================
      ! Single-pass path: whole CSR + lmat scratch fits on the device.
      ! =====================================================================
      !$ACC DATA &
      !$ACC& COPYIN(csr_ptr, csr_id, csr_dist, weights, bintable6, region, &
      !$ACC&        hub_perm, reg_lo, gang_grp, gang_lo, gang_n) &
      !$ACC& CREATE(lmat_g) COPY(part_n4, part_r4, N4jk, R4jk)

      ! One gang per ig slot; hubs dealt round-robin so gang loads even out.
      ! Each gang exclusively owns lmat_g(:,:,ig) and part_*(:,ig).
      ! VECTOR_LENGTH(64): inner-loop trip counts shrink with k2/k3, so shorter
      ! vectors keep lanes busier than 128 while still filling the SMs.
      !$ACC PARALLEL LOOP GANG VECTOR_LENGTH(64) &
      !$ACC& PRIVATE(ig, i, p, g, k1, k2, nn_i, id1, id2, base_i, &
      !$ACC&         ind1, ind2, ind4, w12, rand12, jr1, jr2, jr3, &
      !$ACC&         mixed, dojk, do2, do3)
      do ig = 1, ngang
        ! Gang ig belongs to group g: hubs hub_perm(reg_lo(g):reg_lo(g+1)-1)
        ! dealt round-robin over the gang_n(g) gangs of the group.
        g = gang_grp(ig)
        do p = reg_lo(g) + ig - gang_lo(g), reg_lo(g + 1) - 1, gang_n(g)
          i = hub_perm(p)
          nn_i = int(csr_ptr(i + 1) - csr_ptr(i))
          if (nn_i <= 2) cycle

          base_i = csr_ptr(i) - 1
          jr1 = 0
          mixed = .false.
          if (njk_g > 0) jr1 = region(i)

          ! ---- Phase 1: precompute local connectivity matrix ---------------
          ! C(nn_i, 2) binary searches.  For fixed k1 all lanes search the SAME
          ! id1 adjacency list (different targets) → warp-coherent, L1-cached.
          do k1 = 1, nn_i
            id1 = csr_id(base_i + k1)
            if (njk_g > 0) then
              if (region(id1) /= jr1) mixed = .true.
            end if
            !$ACC LOOP VECTOR PRIVATE(k2, id2)
            do k2 = k1 + 1, nn_i
              id2 = csr_id(base_i + k2)
              lmat_g(k2, k1, ig) = &
                csr_find_dist_bin(csr_ptr, csr_id, csr_dist, id1, id2)
            end do
          end do

          ! Run the jackknife block only for hubs whose neighbours cross a
          ! region boundary (interior hubs are fully covered by the per-gang
          ! partials), or always under the direct-atomic fallback.
          dojk = (njk_g > 0) .and. (mixed .or. direct_hub)

          ! ---- Phase 2: triple loop; k3 across lanes, lmat reads coalesced -
          do k1 = 1, nn_i
            ind1 = csr_dist(base_i + k1)
            id1  = csr_id  (base_i + k1)
            if (dojk) then
              jr2 = region(id1)
              do2 = (jr2 > 0 .and. jr2 /= jr1)
            end if

            do k2 = k1 + 1, nn_i
              ind4 = lmat_g(k2, k1, ig)   ! O(1) lookup: dist bin of k1→k2
              if (ind4 == 0) cycle
              ind2 = csr_dist(base_i + k2)
              id2  = csr_id  (base_i + k2)
              if (dojk) then
                jr3 = region(id2)
                do3 = (jr3 > 0 .and. jr3 /= jr1 .and. jr3 /= jr2)
              end if
              w12  = weights(i) * weights(id1) * weights(id2)
              rand12 = (i > num_data_g .and. id1 > num_data_g .and. &
                        id2 > num_data_g)

              !$ACC LOOP VECTOR PRIVATE(k3, id3, config_idx, ind3, ind5, ind6, w4, jr4, do4)
              do k3 = k2 + 1, nn_i
                ind5 = lmat_g(k3, k1, ig)   ! stride-1 across lanes
                if (ind5 == 0) cycle
                ind6 = lmat_g(k3, k2, ig)   ! stride-1 across lanes
                if (ind6 == 0) cycle

                ind3 = csr_dist(base_i + k3)

                config_idx = abs(bintable6(int(ind1), int(ind2), int(ind3), &
                                           int(ind4), int(ind5), int(ind6)))
                if (config_idx == 0) cycle

                id3 = csr_id(base_i + k3)
                w4  = w12 * weights(id3)

                !$ACC ATOMIC UPDATE
                part_n4(config_idx, ig) = part_n4(config_idx, ig) + w4

                if (rand12 .and. id3 > num_data_g) then
                  !$ACC ATOMIC UPDATE
                  part_r4(config_idx, ig) = part_r4(config_idx, ig) + w4
                end if

                ! Jackknife (see module header / jk_gang_layout): the hub-region
                ! term rides on part_n4/part_r4; only neighbour regions that differ
                ! from the hub's need atomics, and only mixed hubs get here (dojk).
                ! direct_hub: fallback that also does the hub term with atomics.
                if (dojk) then
                  jr4 = region(id3)
                  do4 = (jr4 > 0 .and. jr4 /= jr1 .and. jr4 /= jr2 .and. jr4 /= jr3)
                  if (direct_hub .and. jr1 > 0) then
                    !$ACC ATOMIC UPDATE
                    N4jk(config_idx, 1, jr1) = N4jk(config_idx, 1, jr1) + w4
                  end if
                  if (do2) then
                    !$ACC ATOMIC UPDATE
                    N4jk(config_idx, 1, jr2) = N4jk(config_idx, 1, jr2) + w4
                  end if
                  if (do3) then
                    !$ACC ATOMIC UPDATE
                    N4jk(config_idx, 1, jr3) = N4jk(config_idx, 1, jr3) + w4
                  end if
                  if (do4) then
                    !$ACC ATOMIC UPDATE
                    N4jk(config_idx, 1, jr4) = N4jk(config_idx, 1, jr4) + w4
                  end if
                  if (rand12 .and. id3 > num_data_g) then
                    if (direct_hub .and. jr1 > 0) then
                      !$ACC ATOMIC UPDATE
                      R4jk(config_idx, 1, jr1) = R4jk(config_idx, 1, jr1) + w4
                    end if
                    if (do2) then
                      !$ACC ATOMIC UPDATE
                      R4jk(config_idx, 1, jr2) = R4jk(config_idx, 1, jr2) + w4
                    end if
                    if (do3) then
                      !$ACC ATOMIC UPDATE
                      R4jk(config_idx, 1, jr3) = R4jk(config_idx, 1, jr3) + w4
                    end if
                    if (do4) then
                      !$ACC ATOMIC UPDATE
                      R4jk(config_idx, 1, jr4) = R4jk(config_idx, 1, jr4) + w4
                    end if
                  end if
                end if

              end do  ! k3 (VECTOR)
            end do  ! k2
          end do  ! k1

        end do  ! i  (gang-sequential hubs)
      end do  ! ig  (GANG)
      !$ACC END PARALLEL LOOP
      !$ACC END DATA

    else
      ! =====================================================================
      ! Chunked path: tiles over (hub window) × (search-window pairs cw1<=cw2).
      ! Phase 1a builds lmat columns for neighbors in cw1, Phase 1b for cw2;
      ! Phase 2 processes pairs with id1 in cw1 and id2 in cw2 — each
      ! quadruplet lands in exactly one tile.
      ! =====================================================================
      call csr_make_splits(win_edges, nwin, split_id, split_e)

      maxw = 0
      do cw1 = 1, nwin
        maxw = max(maxw, split_e(cw1 + 1) - split_e(cw1))
      end do
      if (cfg%rank == 0) &
        print '("4PCF GPU chunked: ",i0," windows (",i0," tiles), window <= ",f6.2," GB")', &
              nwin, nwin * nwin * (nwin + 1) / 2, real(maxw, kdkind) * 5.0d0 / 1.0d9

      allocate(stage1_id(maxw), stage1_dist(maxw))
      allocate(stage2_id(maxw), stage2_dist(maxw))

      !$ACC DATA &
      !$ACC& COPYIN(csr_ptr, weights, bintable6, region) &
      !$ACC& CREATE(lmat_g) COPY(part_n4, part_r4, N4jk, R4jk)

      do cw1 = 1, nwin
        w1lo = split_id(cw1)
        w1hi = split_id(cw1 + 1) - 1
        off1 = split_e(cw1) - 1
        n1   = split_e(cw1 + 1) - split_e(cw1)
        stage1_id  (1:n1) = csr_id  (split_e(cw1):split_e(cw1 + 1) - 1)
        stage1_dist(1:n1) = csr_dist(split_e(cw1):split_e(cw1 + 1) - 1)

        !$ACC DATA COPYIN(stage1_id(1:n1), stage1_dist(1:n1))
        do hw = 1, nwin
          hublo = max(istart, split_id(hw))
          hubhi = min(iend,   split_id(hw + 1) - 1)
          if (hublo > hubhi) cycle
          elo_h = split_e(hw)
          ehi_h = split_e(hw + 1) - 1
          ! Gang/region layout for THIS hub window.  Catalogue order tends
          ! to follow sky position, so with a single global layout most
          ! region groups would sit idle on any one window; the window's
          ! partials are folded on the host after its cw2 loop.
          call jk_gang_layout(hublo, hubhi, ngang, njk_g, hub_perm, reg_lo, &
                              gang_grp, gang_lo, gang_n, direct_hub, quiet=.true.)

          !$ACC DATA COPYIN(csr_id(elo_h:ehi_h), csr_dist(elo_h:ehi_h), &
          !$ACC&            hub_perm, reg_lo, gang_grp, gang_lo, gang_n)
          do cw2 = cw1, nwin
            w2lo = split_id(cw2)
            w2hi = split_id(cw2 + 1) - 1
            off2 = split_e(cw2) - 1
            if (cw2 /= cw1) then
              do1b = 1
              n2 = split_e(cw2 + 1) - split_e(cw2)
              stage2_id  (1:n2) = csr_id  (split_e(cw2):split_e(cw2 + 1) - 1)
              stage2_dist(1:n2) = csr_dist(split_e(cw2):split_e(cw2 + 1) - 1)
            else
              do1b = 0
              n2 = 1   ! minimal dummy mapping; kernel never reads it
            end if

            !$ACC DATA COPYIN(stage2_id(1:n2), stage2_dist(1:n2))
            !$ACC PARALLEL LOOP GANG VECTOR_LENGTH(64) &
            !$ACC& PRIVATE(ig, i, p, g, k1, k2, nn_i, id1, id2, base_i, &
            !$ACC&         ind1, ind2, ind4, w12, rand12, jr1, jr2, jr3, &
            !$ACC&         mixed, dojk, do2, do3)
            do ig = 1, ngang
              g = gang_grp(ig)
              do p = reg_lo(g) + ig - gang_lo(g), reg_lo(g + 1) - 1, gang_n(g)
                i = hub_perm(p)
                if (i < hublo .or. i > hubhi) cycle
                nn_i = int(csr_ptr(i + 1) - csr_ptr(i))
                if (nn_i <= 2) cycle

                base_i = csr_ptr(i) - 1
                jr1 = 0
                mixed = .false.
                if (njk_g > 0) jr1 = region(i)

                ! ---- Phase 1a: lmat columns for neighbors in window cw1 ----
                do k1 = 1, nn_i
                  id1 = csr_id(base_i + k1)
                  if (njk_g > 0) then
                    if (region(id1) /= jr1) mixed = .true.
                  end if
                  if (id1 < w1lo .or. id1 > w1hi) cycle
                  !$ACC LOOP VECTOR PRIVATE(k2, id2)
                  do k2 = k1 + 1, nn_i
                    id2 = csr_id(base_i + k2)
                    lmat_g(k2, k1, ig) = csr_find_dist_bin_win( &
                      csr_ptr, stage1_id, stage1_dist, off1, id1, id2)
                  end do
                end do

                ! ---- Phase 1b: lmat columns for neighbors in window cw2 ----
                if (do1b == 1) then
                  do k1 = 1, nn_i
                    id1 = csr_id(base_i + k1)
                    if (id1 < w2lo .or. id1 > w2hi) cycle
                    !$ACC LOOP VECTOR PRIVATE(k2, id2)
                    do k2 = k1 + 1, nn_i
                      id2 = csr_id(base_i + k2)
                      lmat_g(k2, k1, ig) = csr_find_dist_bin_win( &
                        csr_ptr, stage2_id, stage2_dist, off2, id1, id2)
                    end do
                  end do
                end if

                ! Run the jackknife block only for hubs whose neighbours cross a
                ! region boundary (interior hubs are fully covered by the per-gang
                ! partials), or always under the direct-atomic fallback.
                dojk = (njk_g > 0) .and. (mixed .or. direct_hub)

                ! ---- Phase 2: pairs with id1 in cw1, id2 in cw2 ------------
                do k1 = 1, nn_i
                  id1 = csr_id(base_i + k1)
                  if (id1 < w1lo .or. id1 > w1hi) cycle
                  ind1 = csr_dist(base_i + k1)
                  if (dojk) then
                    jr2 = region(id1)
                    do2 = (jr2 > 0 .and. jr2 /= jr1)
                  end if

                  do k2 = k1 + 1, nn_i
                    id2 = csr_id(base_i + k2)
                    if (id2 < w2lo .or. id2 > w2hi) cycle
                    ind4 = lmat_g(k2, k1, ig)
                    if (ind4 == 0) cycle
                    ind2 = csr_dist(base_i + k2)
                    if (dojk) then
                      jr3 = region(id2)
                      do3 = (jr3 > 0 .and. jr3 /= jr1 .and. jr3 /= jr2)
                    end if
                    w12  = weights(i) * weights(id1) * weights(id2)
                    rand12 = (i > num_data_g .and. id1 > num_data_g .and. &
                              id2 > num_data_g)

                    !$ACC LOOP VECTOR PRIVATE(k3, id3, config_idx, ind3, ind5, ind6, w4, jr4, do4)
                    do k3 = k2 + 1, nn_i
                      ind5 = lmat_g(k3, k1, ig)
                      if (ind5 == 0) cycle
                      ind6 = lmat_g(k3, k2, ig)
                      if (ind6 == 0) cycle

                      ind3 = csr_dist(base_i + k3)

                      config_idx = abs(bintable6(int(ind1), int(ind2), int(ind3), &
                                                 int(ind4), int(ind5), int(ind6)))
                      if (config_idx == 0) cycle

                      id3 = csr_id(base_i + k3)
                      w4  = w12 * weights(id3)

                      !$ACC ATOMIC UPDATE
                      part_n4(config_idx, ig) = part_n4(config_idx, ig) + w4

                      if (rand12 .and. id3 > num_data_g) then
                        !$ACC ATOMIC UPDATE
                        part_r4(config_idx, ig) = part_r4(config_idx, ig) + w4
                      end if

                      ! Jackknife (see module header / jk_gang_layout): the hub-region
                      ! term rides on part_n4/part_r4; only neighbour regions that differ
                      ! from the hub's need atomics, and only mixed hubs get here (dojk).
                      ! direct_hub: fallback that also does the hub term with atomics.
                      if (dojk) then
                        jr4 = region(id3)
                        do4 = (jr4 > 0 .and. jr4 /= jr1 .and. jr4 /= jr2 .and. jr4 /= jr3)
                        if (direct_hub .and. jr1 > 0) then
                          !$ACC ATOMIC UPDATE
                          N4jk(config_idx, 1, jr1) = N4jk(config_idx, 1, jr1) + w4
                        end if
                        if (do2) then
                          !$ACC ATOMIC UPDATE
                          N4jk(config_idx, 1, jr2) = N4jk(config_idx, 1, jr2) + w4
                        end if
                        if (do3) then
                          !$ACC ATOMIC UPDATE
                          N4jk(config_idx, 1, jr3) = N4jk(config_idx, 1, jr3) + w4
                        end if
                        if (do4) then
                          !$ACC ATOMIC UPDATE
                          N4jk(config_idx, 1, jr4) = N4jk(config_idx, 1, jr4) + w4
                        end if
                        if (rand12 .and. id3 > num_data_g) then
                          if (direct_hub .and. jr1 > 0) then
                            !$ACC ATOMIC UPDATE
                            R4jk(config_idx, 1, jr1) = R4jk(config_idx, 1, jr1) + w4
                          end if
                          if (do2) then
                            !$ACC ATOMIC UPDATE
                            R4jk(config_idx, 1, jr2) = R4jk(config_idx, 1, jr2) + w4
                          end if
                          if (do3) then
                            !$ACC ATOMIC UPDATE
                            R4jk(config_idx, 1, jr3) = R4jk(config_idx, 1, jr3) + w4
                          end if
                          if (do4) then
                            !$ACC ATOMIC UPDATE
                            R4jk(config_idx, 1, jr4) = R4jk(config_idx, 1, jr4) + w4
                          end if
                        end if
                      end if

                    end do  ! k3 (VECTOR)
                  end do  ! k2
                end do  ! k1

              end do  ! i  (gang-sequential hubs)
            end do  ! ig  (GANG)
            !$ACC END PARALLEL LOOP
            !$ACC END DATA
          end do  ! cw2
          !$ACC END DATA

          ! Fold this hub window's per-gang partials on the host (the
          ! group -> region mapping is per window) and re-zero them on the
          ! device for the next window.
          !$ACC UPDATE HOST(part_n4, part_r4)
          call fold_4pcf_partials(n_configs_g, 1, ngang, part_n4, part_r4, &
                                  gang_grp, njk_g > 0 .and. .not. direct_hub, &
                                  hub_n4jk, hub_r4jk)
          !$ACC PARALLEL LOOP COLLAPSE(2) PRESENT(part_n4, part_r4)
          do ig = 1, ngang
            do config_idx = 1, n_configs_g
              part_n4(config_idx, ig) = 0.0d0
              part_r4(config_idx, ig) = 0.0d0
            end do
          end do
        end do  ! hw
        !$ACC END DATA
      end do  ! cw1
      !$ACC END DATA
      folded = .true.

      deallocate(stage1_id, stage1_dist, stage2_id, stage2_dist)
      deallocate(split_id, split_e)
    end if

    ! Single pass: fold the partials once (the chunked path folded them per
    ! hub window).  Then add the hub-region jackknife terms to N4jk/R4jk,
    ! which now hold the neighbour-region atomics copied back from the device.
    if (.not. folded) &
      call fold_4pcf_partials(n_configs_g, 1, ngang, part_n4, part_r4, &
                              gang_grp, njk_g > 0 .and. .not. direct_hub, &
                              hub_n4jk, hub_r4jk)
    if (njk_g > 0) then
      N4jk(:, 1:1, :) = N4jk(:, 1:1, :) + hub_n4jk
      R4jk(:, 1:1, :) = R4jk(:, 1:1, :) + hub_r4jk
    end if

    deallocate(lmat_g, part_n4, part_r4, hub_n4jk, hub_r4jk)
    deallocate(hub_perm, reg_lo, gang_grp, gang_lo, gang_n)

    call finish_4pcf_output(.false., istart, iend)
  end subroutine query_graph_4pcf_gpu


  ! ---------------------------------------------------------------------------
  ! 4PCF with parity decomposition — same lmat + gang-slot + pair-window-tile
  ! design as query_graph_4pcf_gpu, with two accumulator channels per config
  ! (even = plain w4, odd = parity_flip*sign_V*w4).  The direction pixels
  ! p1/p2/p3 all come from the HUB's own row (csr_phi), so the chunked hub
  ! window carries csr_phi sections while the staged search windows need only
  ! id+dist.  Window sizing uses 4 resident-window budgets to absorb the extra
  ! two bytes/edge of the hub window's int16 phi section (7+5+5 <= 4x5).
  ! ---------------------------------------------------------------------------
  subroutine query_graph_4pcf_parity_gpu(istart, iend)
    integer, intent(in) :: istart, iend

    integer :: ig, i, k1, k2, k3, nn_i, id1, id2, id3
    integer :: config_idx, parity_flip, sign_V
    integer(int64) :: base_i
    integer(int8) :: ind1, ind2, ind3, ind4, ind5, ind6
    integer :: p1, p2, p3, raw_bin
    real(kdkind) :: w12, w4, vol, sv4
    real(kdkind) :: u1x, u1y, u1z, u2x, u2y, u2z, u3x, u3y, u3z, rn1, rn3
    logical :: rand12, exact_g
    integer :: num_data_g, n_configs_g
    ! Jackknife: hub/neighbor region labels (jr1..jr3 hoisted per loop
    ! level and gang-private; jr4 is read in the vector loop).  mixed: some
    ! neighbour of the hub lies in another region; dojk: run the in-kernel
    ! jackknife block for this hub; do2/do3/do4: neighbour region is new.
    integer :: njk_g, jr1, jr2, jr3, jr4
    logical :: mixed, dojk, do2, do3, do4, direct_hub
    ! Gang/region layout (jk_gang_layout): p indexes hub_perm, g the group.
    integer :: p, g
    integer, allocatable :: hub_perm(:), reg_lo(:), gang_grp(:), gang_lo(:), gang_n(:)
    ! Hub-region jackknife terms folded from the per-gang partials on the
    ! host (kept apart from N4jk/R4jk, whose device copies collect the
    ! neighbour-region atomics and are copied back at the end of the DATA
    ! region); folded = partials already folded per hub window (chunked).
    real(kdkind), allocatable :: hub_n4jk(:,:,:), hub_r4jk(:,:,:)
    logical :: folded

    ! Runtime sizing
    integer :: max_nn, ngang, nwin, env_stat
    integer(int64) :: mnn2, lmat_bytes, freeb, win_edges, jk_bytes
    character(32) :: env

    ! Chunking bookkeeping
    integer :: cw1, cw2, hw, do1b
    integer :: w1lo, w1hi, w2lo, w2hi, hublo, hubhi
    integer(int64) :: off1, off2, n1, n2, maxw, elo_h, ehi_h
    integer(int64) :: pelo, pehi   ! csr_phi window (degenerate under -exactparity)
    integer, allocatable :: split_id(:)
    integer(int64), allocatable :: split_e(:)
    integer, allocatable :: stage1_id(:), stage2_id(:)
    integer(int8), allocatable :: stage1_dist(:), stage2_dist(:)

    integer(int8), allocatable :: lmat_g(:,:,:)
    ! Per-gang spoke unit vectors for -exactparity
    real(kdkind), allocatable :: ug(:,:,:)
    ! Per-gang partials, channel 1 = even, channel 2 = parity-odd.
    real(kdkind), allocatable :: part_n4(:,:,:), part_r4(:,:,:)

    if (.not. cfg%half_graph) then
      print *, 'ERROR: GPU 4PCF parity requires half_graph=.true.'
      stop 1
    end if

    if (cfg%rank == 0) print *, 'Performing 4PCF parity (all configs, GPU lmat)'
    if (cfg%rank == 0) print *, 'begin querying the graph'

    num_data_g  = cfg%num_data
    n_configs_g = cfg%n_configs_4pcf
    exact_g     = cfg%exact_parity
    njk_g       = cfg%njk
    ! dummies so the COPYIN clauses are always valid
    if (.not. allocated(px)) allocate(px(1), py(1), pz(1))
    if (.not. allocated(csr_phi)) allocate(csr_phi(1))

    ! ---- Runtime sizing (see query_graph_4pcf_gpu) ----------------------
    max_nn = max(csr_max_row_len(), 1)
    mnn2   = int(max_nn, int64)**2
    freeb  = gpu_free_mem_bytes()

    ! Jackknife accumulators (N4jk + R4jk, two channels) resident on the
    ! device for the whole query; charged against the budgets below.
    jk_bytes = 0_int64
    if (njk_g > 0) then
      jk_bytes = 2_int64 * int(n_configs_g, int64) * 2_int64 * &
                 int(njk_g, int64) * 8_int64
      if (cfg%rank == 0) &
        print '(" 4PCFp GPU jackknife accumulators: ",f8.2," GB on device")', &
          real(jk_bytes, kdkind) / 1.0d9
      if (freeb >= 0 .and. jk_bytes > max(freeb - RESERVE_BYTES, 0_int64) / 2) then
        print *, 'ERROR: 4PCF jackknife accumulators exceed half the free'
        print *, '       device memory; reduce -njk or -nbins'
        stop 1
      end if
    end if

    if (freeb < 0) then
      ngang = 1024
    else
      ngang = int(min(4096_int64, max(512_int64, &
                      ((freeb - RESERVE_BYTES - jk_bytes) / 3) / mnn2)))
    end if
    lmat_bytes = int(ngang, int64) * mnn2

    ! 4 window budgets: hub window costs 7 B/edge (id 4 + dist 1 + phi 2), the
    ! two staged search windows 5 B/edge — a 4×5 B budget covers 7+5+5 with
    ! slack.
    win_edges = csr_edge_window(4, lmat_bytes + jk_bytes)
    nwin = int((csr_total_edges - 1) / win_edges) + 1
    call get_environment_variable('GRAMSCI_GPU_WIN_EDGES', env, status=env_stat)
    ! Collapse to a single pass only if the resident window really fits: 7 B/edge
    ! against the 20 B/edge budget win_edges was derived from.  The previous
    ! `nwin <= 3` test admitted up to 6 B/edge and over-commits at int16 phi.
    if (env_stat /= 0 .and. nwin > 1 .and. &
        csr_total_edges * 7_int64 <= win_edges * 20_int64) nwin = 1

    if (cfg%rank == 0) &
      print '("4PCFp GPU: max_nn=",i0,"  gangs=",i0,"  lmat scratch=",f6.2," GB",'// &
            '"  windows=",i0)', max_nn, ngang, real(lmat_bytes, kdkind)/1.0d9, nwin

    allocate(lmat_g(max_nn, max_nn, ngang))
    if (exact_g) then
      allocate(ug(max_nn, 3, ngang))
    else
      allocate(ug(1, 3, 1))
    end if
    call jk_gang_layout(istart, iend, ngang, njk_g, hub_perm, reg_lo, &
                        gang_grp, gang_lo, gang_n, direct_hub)
    allocate(hub_n4jk(n_configs_g, 2, max(njk_g, 1)))
    allocate(hub_r4jk(n_configs_g, 2, max(njk_g, 1)))
    hub_n4jk = 0.0d0
    hub_r4jk = 0.0d0
    folded = .false.

    allocate(part_n4(n_configs_g, 2, ngang))
    allocate(part_r4(n_configs_g, 2, ngang))
    part_n4 = 0.0d0
    part_r4 = 0.0d0

    if (nwin == 1) then
      ! =====================================================================
      ! Single-pass path
      ! =====================================================================
      !$ACC DATA &
      !$ACC& COPYIN(csr_ptr, csr_id, csr_dist, csr_phi, weights, &
      !$ACC&        bintable6, chiral_4pcf, dir_x, dir_y, dir_z, px, py, pz, region, &
      !$ACC&        hub_perm, reg_lo, gang_grp, gang_lo, gang_n) &
      !$ACC& CREATE(lmat_g, ug) COPY(part_n4, part_r4, N4jk, R4jk)

      !$ACC PARALLEL LOOP GANG VECTOR_LENGTH(64) &
      !$ACC& PRIVATE(ig, i, p, g, k1, k2, nn_i, id1, id2, base_i, &
      !$ACC&         ind1, ind2, ind4, p1, p2, w12, rand12, &
      !$ACC&         u1x, u1y, u1z, u2x, u2y, u2z, rn1, jr1, jr2, jr3, &
      !$ACC&         mixed, dojk, do2, do3)
      do ig = 1, ngang
        g = gang_grp(ig)
        do p = reg_lo(g) + ig - gang_lo(g), reg_lo(g + 1) - 1, gang_n(g)
          i = hub_perm(p)
          nn_i = int(csr_ptr(i + 1) - csr_ptr(i))
          if (nn_i <= 2) cycle

          base_i = csr_ptr(i) - 1
          jr1 = 0
          mixed = .false.
          if (njk_g > 0) jr1 = region(i)

          ! ---- Phase 1: connectivity matrix ------------------------------
          do k1 = 1, nn_i
            id1 = csr_id(base_i + k1)
            if (njk_g > 0) then
              if (region(id1) /= jr1) mixed = .true.
            end if
            !$ACC LOOP VECTOR PRIVATE(k2, id2)
            do k2 = k1 + 1, nn_i
              id2 = csr_id(base_i + k2)
              lmat_g(k2, k1, ig) = &
                csr_find_dist_bin(csr_ptr, csr_id, csr_dist, id1, id2)
            end do
          end do

          ! ---- Phase 1p: spoke unit vectors for -exactparity -------------
          if (exact_g) then
            !$ACC LOOP VECTOR PRIVATE(k1, id1, u3x, u3y, u3z, rn3)
            do k1 = 1, nn_i
              id1 = csr_id(base_i + k1)
              u3x = px(id1) - px(i)
              u3y = py(id1) - py(i)
              u3z = pz(id1) - pz(i)
              rn3 = 1.0d0 / sqrt(u3x*u3x + u3y*u3y + u3z*u3z)
              ug(k1, 1, ig) = u3x * rn3
              ug(k1, 2, ig) = u3y * rn3
              ug(k1, 3, ig) = u3z * rn3
            end do
          end if

          ! Run the jackknife block only for hubs whose neighbours cross a
          ! region boundary (interior hubs are fully covered by the per-gang
          ! partials), or always under the direct-atomic fallback.
          dojk = (njk_g > 0) .and. (mixed .or. direct_hub)

          ! ---- Phase 2: triple loop with parity channel ------------------
          do k1 = 1, nn_i
            ind1 = csr_dist(base_i + k1)
            id1  = csr_id  (base_i + k1)
            if (dojk) then
              jr2 = region(id1)
              do2 = (jr2 > 0 .and. jr2 /= jr1)
            end if
            if (exact_g) then
              u1x = ug(k1, 1, ig); u1y = ug(k1, 2, ig); u1z = ug(k1, 3, ig)
            else
              p1 = int(csr_phi(base_i + k1))
            end if

            do k2 = k1 + 1, nn_i
              ind4 = lmat_g(k2, k1, ig)
              if (ind4 == 0) cycle
              ind2 = csr_dist(base_i + k2)
              id2  = csr_id  (base_i + k2)
              if (dojk) then
                jr3 = region(id2)
                do3 = (jr3 > 0 .and. jr3 /= jr1 .and. jr3 /= jr2)
              end if
              if (exact_g) then
                u2x = ug(k2, 1, ig); u2y = ug(k2, 2, ig); u2z = ug(k2, 3, ig)
              else
                p2 = int(csr_phi(base_i + k2))
              end if
              w12  = weights(i) * weights(id1) * weights(id2)
              rand12 = (i > num_data_g .and. id1 > num_data_g .and. &
                        id2 > num_data_g)

              !$ACC LOOP VECTOR PRIVATE(k3, id3, config_idx, raw_bin, &
              !$ACC&   parity_flip, sign_V, ind3, ind5, ind6, p3, w4, vol, &
              !$ACC&   u3x, u3y, u3z, rn3, sv4, jr4, do4)
              do k3 = k2 + 1, nn_i
                ind5 = lmat_g(k3, k1, ig)
                if (ind5 == 0) cycle
                ind6 = lmat_g(k3, k2, ig)
                if (ind6 == 0) cycle

                ind3 = csr_dist(base_i + k3)

                raw_bin    = bintable6(int(ind1), int(ind2), int(ind3), &
                                       int(ind4), int(ind5), int(ind6))
                config_idx = abs(raw_bin)
                if (config_idx == 0) cycle
                parity_flip = sign(1, raw_bin) * int(chiral_4pcf(config_idx))

                id3 = csr_id(base_i + k3)
                if (exact_g) then
                  u3x = ug(k3, 1, ig); u3y = ug(k3, 2, ig); u3z = ug(k3, 3, ig)
                  vol = u1x * (u2y*u3z - u2z*u3y) &
                      + u1y * (u2z*u3x - u2x*u3z) &
                      + u1z * (u2x*u3y - u2y*u3x)
                else
                  p3 = int(csr_phi(base_i + k3))
                  vol = dir_x(p1) * (dir_y(p2)*dir_z(p3) - dir_z(p2)*dir_y(p3)) &
                      + dir_y(p1) * (dir_z(p2)*dir_x(p3) - dir_x(p2)*dir_z(p3)) &
                      + dir_z(p1) * (dir_x(p2)*dir_y(p3) - dir_y(p2)*dir_x(p3))
                end if
                if (abs(vol) < VOL_DEGEN_TOL) then
                  sign_V = 0   ! degenerate: no chirality, odd channel gets 0
                else if (vol > 0.0d0) then
                  sign_V = 1
                else
                  sign_V = -1
                end if

                w4  = w12 * weights(id3)
                sv4 = (parity_flip * sign_V) * w4

                !$ACC ATOMIC UPDATE
                part_n4(config_idx, 1, ig) = part_n4(config_idx, 1, ig) + w4
                !$ACC ATOMIC UPDATE
                part_n4(config_idx, 2, ig) = part_n4(config_idx, 2, ig) + sv4

                if (rand12 .and. id3 > num_data_g) then
                  !$ACC ATOMIC UPDATE
                  part_r4(config_idx, 1, ig) = part_r4(config_idx, 1, ig) + w4
                  !$ACC ATOMIC UPDATE
                  part_r4(config_idx, 2, ig) = part_r4(config_idx, 2, ig) + sv4
                end if

                ! Jackknife (see module header / jk_gang_layout): the hub-region
                ! term rides on part_n4/part_r4; only neighbour regions that differ
                ! from the hub's need atomics, and only mixed hubs get here (dojk).
                ! direct_hub: fallback that also does the hub term with atomics.
                if (dojk) then
                  jr4 = region(id3)
                  do4 = (jr4 > 0 .and. jr4 /= jr1 .and. jr4 /= jr2 .and. jr4 /= jr3)
                  if (direct_hub .and. jr1 > 0) then
                    !$ACC ATOMIC UPDATE
                    N4jk(config_idx, 1, jr1) = N4jk(config_idx, 1, jr1) + w4
                    !$ACC ATOMIC UPDATE
                    N4jk(config_idx, 2, jr1) = N4jk(config_idx, 2, jr1) + sv4
                  end if
                  if (do2) then
                    !$ACC ATOMIC UPDATE
                    N4jk(config_idx, 1, jr2) = N4jk(config_idx, 1, jr2) + w4
                    !$ACC ATOMIC UPDATE
                    N4jk(config_idx, 2, jr2) = N4jk(config_idx, 2, jr2) + sv4
                  end if
                  if (do3) then
                    !$ACC ATOMIC UPDATE
                    N4jk(config_idx, 1, jr3) = N4jk(config_idx, 1, jr3) + w4
                    !$ACC ATOMIC UPDATE
                    N4jk(config_idx, 2, jr3) = N4jk(config_idx, 2, jr3) + sv4
                  end if
                  if (do4) then
                    !$ACC ATOMIC UPDATE
                    N4jk(config_idx, 1, jr4) = N4jk(config_idx, 1, jr4) + w4
                    !$ACC ATOMIC UPDATE
                    N4jk(config_idx, 2, jr4) = N4jk(config_idx, 2, jr4) + sv4
                  end if
                  if (rand12 .and. id3 > num_data_g) then
                    if (direct_hub .and. jr1 > 0) then
                      !$ACC ATOMIC UPDATE
                      R4jk(config_idx, 1, jr1) = R4jk(config_idx, 1, jr1) + w4
                      !$ACC ATOMIC UPDATE
                      R4jk(config_idx, 2, jr1) = R4jk(config_idx, 2, jr1) + sv4
                    end if
                    if (do2) then
                      !$ACC ATOMIC UPDATE
                      R4jk(config_idx, 1, jr2) = R4jk(config_idx, 1, jr2) + w4
                      !$ACC ATOMIC UPDATE
                      R4jk(config_idx, 2, jr2) = R4jk(config_idx, 2, jr2) + sv4
                    end if
                    if (do3) then
                      !$ACC ATOMIC UPDATE
                      R4jk(config_idx, 1, jr3) = R4jk(config_idx, 1, jr3) + w4
                      !$ACC ATOMIC UPDATE
                      R4jk(config_idx, 2, jr3) = R4jk(config_idx, 2, jr3) + sv4
                    end if
                    if (do4) then
                      !$ACC ATOMIC UPDATE
                      R4jk(config_idx, 1, jr4) = R4jk(config_idx, 1, jr4) + w4
                      !$ACC ATOMIC UPDATE
                      R4jk(config_idx, 2, jr4) = R4jk(config_idx, 2, jr4) + sv4
                    end if
                  end if
                end if

              end do  ! k3 (VECTOR)
            end do  ! k2
          end do  ! k1

        end do  ! i
      end do  ! ig (GANG)
      !$ACC END PARALLEL LOOP
      !$ACC END DATA

    else
      ! =====================================================================
      ! Chunked path: (hub window) × (search-window pairs cw1<=cw2) tiles.
      ! =====================================================================
      call csr_make_splits(win_edges, nwin, split_id, split_e)

      maxw = 0
      do cw1 = 1, nwin
        maxw = max(maxw, split_e(cw1 + 1) - split_e(cw1))
      end do
      if (cfg%rank == 0) &
        print '("4PCFp GPU chunked: ",i0," windows (",i0," tiles), window <= ",f6.2," GB")', &
              nwin, nwin * nwin * (nwin + 1) / 2, real(maxw, kdkind) * 7.0d0 / 1.0d9

      allocate(stage1_id(maxw), stage1_dist(maxw))
      allocate(stage2_id(maxw), stage2_dist(maxw))

      !$ACC DATA &
      !$ACC& COPYIN(csr_ptr, weights, bintable6, chiral_4pcf, dir_x, dir_y, dir_z, &
      !$ACC&        px, py, pz, region) &
      !$ACC& CREATE(lmat_g, ug) COPY(part_n4, part_r4, N4jk, R4jk)

      do cw1 = 1, nwin
        w1lo = split_id(cw1)
        w1hi = split_id(cw1 + 1) - 1
        off1 = split_e(cw1) - 1
        n1   = split_e(cw1 + 1) - split_e(cw1)
        stage1_id  (1:n1) = csr_id  (split_e(cw1):split_e(cw1 + 1) - 1)
        stage1_dist(1:n1) = csr_dist(split_e(cw1):split_e(cw1 + 1) - 1)

        !$ACC DATA COPYIN(stage1_id(1:n1), stage1_dist(1:n1))
        do hw = 1, nwin
          hublo = max(istart, split_id(hw))
          hubhi = min(iend,   split_id(hw + 1) - 1)
          if (hublo > hubhi) cycle
          elo_h = split_e(hw)
          ehi_h = split_e(hw + 1) - 1
          ! Gang/region layout for THIS hub window.  Catalogue order tends
          ! to follow sky position, so with a single global layout most
          ! region groups would sit idle on any one window; the window's
          ! partials are folded on the host after its cw2 loop.
          call jk_gang_layout(hublo, hubhi, ngang, njk_g, hub_perm, reg_lo, &
                              gang_grp, gang_lo, gang_n, direct_hub, quiet=.true.)
          ! csr_phi is a 1-element dummy under -exactparity (not stored)
          if (exact_g) then
            pelo = 1_int64
            pehi = 1_int64
          else
            pelo = elo_h
            pehi = ehi_h
          end if

          !$ACC DATA COPYIN(csr_id(elo_h:ehi_h), csr_dist(elo_h:ehi_h), &
          !$ACC&            csr_phi(pelo:pehi), &
          !$ACC&            hub_perm, reg_lo, gang_grp, gang_lo, gang_n)
          do cw2 = cw1, nwin
            w2lo = split_id(cw2)
            w2hi = split_id(cw2 + 1) - 1
            off2 = split_e(cw2) - 1
            if (cw2 /= cw1) then
              do1b = 1
              n2 = split_e(cw2 + 1) - split_e(cw2)
              stage2_id  (1:n2) = csr_id  (split_e(cw2):split_e(cw2 + 1) - 1)
              stage2_dist(1:n2) = csr_dist(split_e(cw2):split_e(cw2 + 1) - 1)
            else
              do1b = 0
              n2 = 1
            end if

            !$ACC DATA COPYIN(stage2_id(1:n2), stage2_dist(1:n2))
            !$ACC PARALLEL LOOP GANG VECTOR_LENGTH(64) &
            !$ACC& PRIVATE(ig, i, p, g, k1, k2, nn_i, id1, id2, base_i, &
            !$ACC&         ind1, ind2, ind4, p1, p2, w12, rand12, &
            !$ACC&         u1x, u1y, u1z, u2x, u2y, u2z, rn1, jr1, jr2, jr3, &
            !$ACC&         mixed, dojk, do2, do3)
            do ig = 1, ngang
              g = gang_grp(ig)
              do p = reg_lo(g) + ig - gang_lo(g), reg_lo(g + 1) - 1, gang_n(g)
                i = hub_perm(p)
                if (i < hublo .or. i > hubhi) cycle
                nn_i = int(csr_ptr(i + 1) - csr_ptr(i))
                if (nn_i <= 2) cycle

                base_i = csr_ptr(i) - 1
                jr1 = 0
                mixed = .false.
                if (njk_g > 0) jr1 = region(i)

                ! ---- Phase 1a: columns for neighbors in window cw1 ---------
                do k1 = 1, nn_i
                  id1 = csr_id(base_i + k1)
                  if (njk_g > 0) then
                    if (region(id1) /= jr1) mixed = .true.
                  end if
                  if (id1 < w1lo .or. id1 > w1hi) cycle
                  !$ACC LOOP VECTOR PRIVATE(k2, id2)
                  do k2 = k1 + 1, nn_i
                    id2 = csr_id(base_i + k2)
                    lmat_g(k2, k1, ig) = csr_find_dist_bin_win( &
                      csr_ptr, stage1_id, stage1_dist, off1, id1, id2)
                  end do
                end do

                ! ---- Phase 1b: columns for neighbors in window cw2 ---------
                if (do1b == 1) then
                  do k1 = 1, nn_i
                    id1 = csr_id(base_i + k1)
                    if (id1 < w2lo .or. id1 > w2hi) cycle
                    !$ACC LOOP VECTOR PRIVATE(k2, id2)
                    do k2 = k1 + 1, nn_i
                      id2 = csr_id(base_i + k2)
                      lmat_g(k2, k1, ig) = csr_find_dist_bin_win( &
                        csr_ptr, stage2_id, stage2_dist, off2, id1, id2)
                    end do
                  end do
                end if

                ! ---- Phase 1p: spoke unit vectors for -exactparity ---------
                if (exact_g) then
                  !$ACC LOOP VECTOR PRIVATE(k1, id1, u3x, u3y, u3z, rn3)
                  do k1 = 1, nn_i
                    id1 = csr_id(base_i + k1)
                    u3x = px(id1) - px(i)
                    u3y = py(id1) - py(i)
                    u3z = pz(id1) - pz(i)
                    rn3 = 1.0d0 / sqrt(u3x*u3x + u3y*u3y + u3z*u3z)
                    ug(k1, 1, ig) = u3x * rn3
                    ug(k1, 2, ig) = u3y * rn3
                    ug(k1, 3, ig) = u3z * rn3
                  end do
                end if

                ! Run the jackknife block only for hubs whose neighbours cross a
                ! region boundary (interior hubs are fully covered by the per-gang
                ! partials), or always under the direct-atomic fallback.
                dojk = (njk_g > 0) .and. (mixed .or. direct_hub)

                ! ---- Phase 2: pairs with id1 in cw1, id2 in cw2 ------------
                do k1 = 1, nn_i
                  id1 = csr_id(base_i + k1)
                  if (id1 < w1lo .or. id1 > w1hi) cycle
                  ind1 = csr_dist(base_i + k1)
                  if (dojk) then
                    jr2 = region(id1)
                    do2 = (jr2 > 0 .and. jr2 /= jr1)
                  end if
                  if (exact_g) then
                    u1x = ug(k1, 1, ig); u1y = ug(k1, 2, ig); u1z = ug(k1, 3, ig)
                  else
                    p1 = int(csr_phi(base_i + k1))
                  end if

                  do k2 = k1 + 1, nn_i
                    id2 = csr_id(base_i + k2)
                    if (id2 < w2lo .or. id2 > w2hi) cycle
                    ind4 = lmat_g(k2, k1, ig)
                    if (ind4 == 0) cycle
                    ind2 = csr_dist(base_i + k2)
                    if (dojk) then
                      jr3 = region(id2)
                      do3 = (jr3 > 0 .and. jr3 /= jr1 .and. jr3 /= jr2)
                    end if
                    if (exact_g) then
                      u2x = ug(k2, 1, ig); u2y = ug(k2, 2, ig); u2z = ug(k2, 3, ig)
                    else
                      p2 = int(csr_phi(base_i + k2))
                    end if
                    w12  = weights(i) * weights(id1) * weights(id2)
                    rand12 = (i > num_data_g .and. id1 > num_data_g .and. &
                              id2 > num_data_g)

                    !$ACC LOOP VECTOR PRIVATE(k3, id3, config_idx, raw_bin, &
                    !$ACC&   parity_flip, sign_V, ind3, ind5, ind6, p3, w4, vol, &
                    !$ACC&   u3x, u3y, u3z, rn3, sv4, jr4, do4)
                    do k3 = k2 + 1, nn_i
                      ind5 = lmat_g(k3, k1, ig)
                      if (ind5 == 0) cycle
                      ind6 = lmat_g(k3, k2, ig)
                      if (ind6 == 0) cycle

                      ind3 = csr_dist(base_i + k3)

                      raw_bin    = bintable6(int(ind1), int(ind2), int(ind3), &
                                             int(ind4), int(ind5), int(ind6))
                      config_idx = abs(raw_bin)
                      if (config_idx == 0) cycle
                      parity_flip = sign(1, raw_bin) * int(chiral_4pcf(config_idx))

                      id3 = csr_id(base_i + k3)
                      if (exact_g) then
                        u3x = ug(k3, 1, ig); u3y = ug(k3, 2, ig); u3z = ug(k3, 3, ig)
                        vol = u1x * (u2y*u3z - u2z*u3y) &
                            + u1y * (u2z*u3x - u2x*u3z) &
                            + u1z * (u2x*u3y - u2y*u3x)
                      else
                        p3 = int(csr_phi(base_i + k3))
                        vol = dir_x(p1) * (dir_y(p2)*dir_z(p3) - dir_z(p2)*dir_y(p3)) &
                            + dir_y(p1) * (dir_z(p2)*dir_x(p3) - dir_x(p2)*dir_z(p3)) &
                            + dir_z(p1) * (dir_x(p2)*dir_y(p3) - dir_y(p2)*dir_x(p3))
                      end if
                      if (abs(vol) < VOL_DEGEN_TOL) then
                        sign_V = 0   ! degenerate: no chirality, odd channel gets 0
                      else if (vol > 0.0d0) then
                        sign_V = 1
                      else
                        sign_V = -1
                      end if

                      w4  = w12 * weights(id3)
                      sv4 = (parity_flip * sign_V) * w4

                      !$ACC ATOMIC UPDATE
                      part_n4(config_idx, 1, ig) = part_n4(config_idx, 1, ig) + w4
                      !$ACC ATOMIC UPDATE
                      part_n4(config_idx, 2, ig) = part_n4(config_idx, 2, ig) + sv4

                      if (rand12 .and. id3 > num_data_g) then
                        !$ACC ATOMIC UPDATE
                        part_r4(config_idx, 1, ig) = part_r4(config_idx, 1, ig) + w4
                        !$ACC ATOMIC UPDATE
                        part_r4(config_idx, 2, ig) = part_r4(config_idx, 2, ig) + sv4
                      end if

                      ! Jackknife (see module header / jk_gang_layout): the hub-region
                      ! term rides on part_n4/part_r4; only neighbour regions that differ
                      ! from the hub's need atomics, and only mixed hubs get here (dojk).
                      ! direct_hub: fallback that also does the hub term with atomics.
                      if (dojk) then
                        jr4 = region(id3)
                        do4 = (jr4 > 0 .and. jr4 /= jr1 .and. jr4 /= jr2 .and. jr4 /= jr3)
                        if (direct_hub .and. jr1 > 0) then
                          !$ACC ATOMIC UPDATE
                          N4jk(config_idx, 1, jr1) = N4jk(config_idx, 1, jr1) + w4
                          !$ACC ATOMIC UPDATE
                          N4jk(config_idx, 2, jr1) = N4jk(config_idx, 2, jr1) + sv4
                        end if
                        if (do2) then
                          !$ACC ATOMIC UPDATE
                          N4jk(config_idx, 1, jr2) = N4jk(config_idx, 1, jr2) + w4
                          !$ACC ATOMIC UPDATE
                          N4jk(config_idx, 2, jr2) = N4jk(config_idx, 2, jr2) + sv4
                        end if
                        if (do3) then
                          !$ACC ATOMIC UPDATE
                          N4jk(config_idx, 1, jr3) = N4jk(config_idx, 1, jr3) + w4
                          !$ACC ATOMIC UPDATE
                          N4jk(config_idx, 2, jr3) = N4jk(config_idx, 2, jr3) + sv4
                        end if
                        if (do4) then
                          !$ACC ATOMIC UPDATE
                          N4jk(config_idx, 1, jr4) = N4jk(config_idx, 1, jr4) + w4
                          !$ACC ATOMIC UPDATE
                          N4jk(config_idx, 2, jr4) = N4jk(config_idx, 2, jr4) + sv4
                        end if
                        if (rand12 .and. id3 > num_data_g) then
                          if (direct_hub .and. jr1 > 0) then
                            !$ACC ATOMIC UPDATE
                            R4jk(config_idx, 1, jr1) = R4jk(config_idx, 1, jr1) + w4
                            !$ACC ATOMIC UPDATE
                            R4jk(config_idx, 2, jr1) = R4jk(config_idx, 2, jr1) + sv4
                          end if
                          if (do2) then
                            !$ACC ATOMIC UPDATE
                            R4jk(config_idx, 1, jr2) = R4jk(config_idx, 1, jr2) + w4
                            !$ACC ATOMIC UPDATE
                            R4jk(config_idx, 2, jr2) = R4jk(config_idx, 2, jr2) + sv4
                          end if
                          if (do3) then
                            !$ACC ATOMIC UPDATE
                            R4jk(config_idx, 1, jr3) = R4jk(config_idx, 1, jr3) + w4
                            !$ACC ATOMIC UPDATE
                            R4jk(config_idx, 2, jr3) = R4jk(config_idx, 2, jr3) + sv4
                          end if
                          if (do4) then
                            !$ACC ATOMIC UPDATE
                            R4jk(config_idx, 1, jr4) = R4jk(config_idx, 1, jr4) + w4
                            !$ACC ATOMIC UPDATE
                            R4jk(config_idx, 2, jr4) = R4jk(config_idx, 2, jr4) + sv4
                          end if
                        end if
                      end if

                    end do  ! k3 (VECTOR)
                  end do  ! k2
                end do  ! k1

              end do  ! i
            end do  ! ig (GANG)
            !$ACC END PARALLEL LOOP
            !$ACC END DATA
          end do  ! cw2
          !$ACC END DATA

          ! Fold this hub window's per-gang partials on the host (the
          ! group -> region mapping is per window) and re-zero them on the
          ! device for the next window.
          !$ACC UPDATE HOST(part_n4, part_r4)
          call fold_4pcf_partials(n_configs_g, 2, ngang, part_n4, part_r4, &
                                  gang_grp, njk_g > 0 .and. .not. direct_hub, &
                                  hub_n4jk, hub_r4jk)
          !$ACC PARALLEL LOOP COLLAPSE(2) PRESENT(part_n4, part_r4)
          do ig = 1, ngang
            do config_idx = 1, n_configs_g
              part_n4(config_idx, 1, ig) = 0.0d0
              part_n4(config_idx, 2, ig) = 0.0d0
              part_r4(config_idx, 1, ig) = 0.0d0
              part_r4(config_idx, 2, ig) = 0.0d0
            end do
          end do
        end do  ! hw
        !$ACC END DATA
      end do  ! cw1
      !$ACC END DATA
      folded = .true.

      deallocate(stage1_id, stage1_dist, stage2_id, stage2_dist)
      deallocate(split_id, split_e)
    end if

    ! Single pass: fold the partials once (the chunked path folded them per
    ! hub window).  Then add the hub-region jackknife terms to N4jk/R4jk,
    ! which now hold the neighbour-region atomics copied back from the device.
    if (.not. folded) &
      call fold_4pcf_partials(n_configs_g, 2, ngang, part_n4, part_r4, &
                              gang_grp, njk_g > 0 .and. .not. direct_hub, &
                              hub_n4jk, hub_r4jk)
    if (njk_g > 0) then
      N4jk(:, 1:2, :) = N4jk(:, 1:2, :) + hub_n4jk
      R4jk(:, 1:2, :) = R4jk(:, 1:2, :) + hub_r4jk
    end if

    deallocate(lmat_g, part_n4, part_r4, hub_n4jk, hub_r4jk)
    deallocate(hub_perm, reg_lo, gang_grp, gang_lo, gang_n)
    if (allocated(ug)) deallocate(ug)

    call finish_4pcf_output(.true., istart, iend)
  end subroutine query_graph_4pcf_parity_gpu

  ! ---------------------------------------------------------------------------
  ! Jackknife gang layout (both 4PCF kernels).
  !
  ! Hubs istart..iend are permuted so that all hubs of one region are
  ! contiguous: hub_perm(reg_lo(g) : reg_lo(g+1)-1) are the hubs of group g,
  ! where group g = r+1 holds region r (g = 1 collects unlabelled hubs,
  ! region 0).  The ngang gangs are split into one contiguous block per
  ! non-empty group, sized in proportion to the group's hub count
  ! (gang_grp(ig) = group of gang ig, gang_lo(g) = first gang of the group,
  ! gang_n(g) = its gang count; hubs are dealt round-robin over the block).
  !
  ! A gang therefore only ever processes hubs of a single region, so its
  ! per-gang partials ARE that region's hub-term contribution to N4jk/R4jk
  ! and are added on the host — no device atomic for the hub region at all.
  ! The jackknife block in the kernel only handles the three neighbours'
  ! regions when they differ from the hub's.
  !
  ! direct_hub = .true. selects the fallback when the non-empty groups
  ! cannot be spread over the gangs with a sensible minimum block size
  ! (more than ngang/4 regions): a single group with every gang, and the
  ! kernel adds the hub term with atomics as well.  Without -njk the same
  ! single-group layout with the identity permutation reproduces the plain
  ! round-robin hub dealing exactly.
  !
  ! The chunked kernels call this once per hub window (quiet=.true.): the
  ! catalogue order usually follows sky position, so one global layout
  ! would leave most region groups idle on any single window.
  ! ---------------------------------------------------------------------------
  subroutine jk_gang_layout(istart, iend, ngang, njk, hub_perm, reg_lo, &
                            gang_grp, gang_lo, gang_n, direct_hub, quiet)
    integer, intent(in) :: istart, iend, ngang, njk
    integer, allocatable, intent(out) :: hub_perm(:), reg_lo(:)
    integer, allocatable, intent(out) :: gang_grp(:), gang_lo(:), gang_n(:)
    logical, intent(out) :: direct_hub
    logical, intent(in), optional :: quiet   ! no diagnostics (per-window calls)

    integer :: nhub, ngrp, nonempty, g, i, ig, r, rem, gbest
    integer(int64) :: nmixed, e
    integer, allocatable :: nh(:), fill(:)
    real(kdkind) :: score, best
    logical :: verbose

    verbose = .true.
    if (present(quiet)) verbose = .not. quiet

    nhub = max(iend - istart + 1, 0)
    ngrp = njk + 1                       ! group g = region g-1
    if (njk <= 0) ngrp = 1

    allocate(nh(ngrp), reg_lo(ngrp + 1), hub_perm(max(nhub, 1)))
    allocate(gang_grp(ngang), gang_lo(ngrp), gang_n(ngrp))

    ! ---- hub counts per group, and the region-sorted permutation --------
    nh = 0
    if (njk > 0) then
      do i = istart, iend
        r = region(i)
        if (r < 0 .or. r > njk) r = 0
        nh(r + 1) = nh(r + 1) + 1
      end do
    else
      nh(1) = nhub
    end if
    nonempty = count(nh > 0)

    ! direct_hub: too many regions to give every non-empty one a block of
    ! gangs (block < 4 gangs) — keep one group and do the hub term with
    ! atomics in the kernel.
    direct_hub = (njk > 0 .and. nonempty * 4 > ngang)
    if (direct_hub .or. njk <= 0) then
      ngrp = 1
      deallocate(nh, reg_lo, gang_lo, gang_n)
      allocate(nh(1), reg_lo(2), gang_lo(1), gang_n(1))
      nh(1) = nhub
    end if

    reg_lo(1) = 1
    do g = 1, ngrp
      reg_lo(g + 1) = reg_lo(g) + nh(g)
    end do
    allocate(fill(ngrp))
    fill = reg_lo(1:ngrp)
    if (ngrp == 1) then
      do i = istart, iend
        hub_perm(i - istart + 1) = i
      end do
    else
      do i = istart, iend                 ! stable counting sort by region
        r = region(i)
        if (r < 0 .or. r > njk) r = 0
        hub_perm(fill(r + 1)) = i
        fill(r + 1) = fill(r + 1) + 1
      end do
    end if

    ! ---- gangs per group, proportional to hub count, >= 1 if non-empty --
    gang_n = 0
    if (ngrp == 1) then
      gang_n(1) = ngang
    else
      do g = 1, ngrp
        if (nh(g) > 0) gang_n(g) = max(1, int((int(ngang, int64) * nh(g)) / nhub))
      end do
      rem = ngang - sum(gang_n)
      do while (rem /= 0)
        ! rem > 0: give a gang to the group with the most hubs per gang;
        ! rem < 0: take one from the group that stays best-off after losing it.
        gbest = 0
        best  = -1.0d0
        do g = 1, ngrp
          if (nh(g) == 0) cycle
          if (rem > 0) then
            score = real(nh(g), kdkind) / real(gang_n(g), kdkind)
          else
            if (gang_n(g) <= 1) cycle
            score = -real(nh(g), kdkind) / real(gang_n(g) - 1, kdkind)
          end if
          if (score > best) then
            best  = score
            gbest = g
          end if
        end do
        if (gbest == 0) exit              ! cannot shrink further (all at 1)
        if (rem > 0) then
          gang_n(gbest) = gang_n(gbest) + 1
          rem = rem - 1
        else
          gang_n(gbest) = gang_n(gbest) - 1
          rem = rem + 1
        end if
      end do
    end if
    ! The split above always distributes exactly ngang gangs: nonempty <=
    ! ngang/4 in this branch, so the shrink loop can always find a group
    ! with more than one gang.  Guard it anyway — an unassigned gang would
    ! otherwise re-walk another gang's hubs and double count.
    if (sum(gang_n) /= ngang) then
      print *, 'ERROR: jk_gang_layout: gang split mismatch ', sum(gang_n), ngang
      stop 1
    end if
    gang_lo(1) = 1
    do g = 1, ngrp - 1
      gang_lo(g + 1) = gang_lo(g) + gang_n(g)
    end do
    do g = 1, ngrp
      do ig = gang_lo(g), gang_lo(g) + gang_n(g) - 1
        gang_grp(ig) = g
      end do
    end do

    if (njk > 0 .and. cfg%rank == 0 .and. verbose) then
      ! Diagnostic: fraction of hubs whose adjacency list crosses a region
      ! boundary — only these run the in-kernel jackknife block.
      nmixed = 0
      !$OMP PARALLEL DO PRIVATE(i, e, r) REDUCTION(+:nmixed) SCHEDULE(STATIC)
      do i = istart, iend
        r = region(i)
        do e = csr_ptr(i), csr_ptr(i + 1) - 1
          if (region(csr_id(e)) /= r) then
            nmixed = nmixed + 1
            exit
          end if
        end do
      end do
      !$OMP END PARALLEL DO
      if (direct_hub) then
        print '(" 4PCF GPU jackknife: ",i0," non-empty regions (of ",i0,") > gangs/4: ",'// &
              '"direct hub-term atomics")', nonempty, njk
      else
        print '(" 4PCF GPU jackknife: ",i0," region gang groups, ",i0,"-",i0," gangs each")', &
          nonempty, minval(gang_n, mask=gang_n > 0), maxval(gang_n)
      end if
      print '("   mixed hubs (neighbours in another region): ",i0," of ",i0," (",f5.1,"%)")', &
        nmixed, nhub, 100.0d0 * real(nmixed, kdkind) / real(max(nhub, 1), kdkind)
    end if

    deallocate(nh, fill)
  end subroutine jk_gang_layout

  ! ---------------------------------------------------------------------------
  ! Fold per-gang partials part_*(n_configs, nch, ngang) into the totals
  ! N4/R4 and, when the gang groups are in use, the hub-region jackknife
  ! accumulators hub_*(:, :, region) (group g of gang_grp holds region g-1;
  ! region 0 = unlabelled hubs contribute to the totals only).  The partials
  ! are passed by sequence association, so the 2-D (n_configs, ngang) arrays
  ! of the non-parity kernel are handled with nch = 1.
  ! ---------------------------------------------------------------------------
  subroutine fold_4pcf_partials(nconf, nch, ngang, part_n4, part_r4, gang_grp, &
                                use_groups, hub_n4jk, hub_r4jk)
    integer, intent(in) :: nconf, nch, ngang, gang_grp(ngang)
    real(kdkind), intent(in) :: part_n4(nconf, nch, ngang), part_r4(nconf, nch, ngang)
    logical, intent(in) :: use_groups
    real(kdkind), intent(inout) :: hub_n4jk(:,:,:), hub_r4jk(:,:,:)
    integer :: c, ch, ig, g

    do ch = 1, nch
      do c = 1, nconf
        N4(c, ch) = N4(c, ch) + sum(part_n4(c, ch, :))
        R4(c, ch) = R4(c, ch) + sum(part_r4(c, ch, :))
      end do
    end do
    if (.not. use_groups) return
    do ig = 1, ngang
      g = gang_grp(ig) - 1
      if (g <= 0) cycle
      hub_n4jk(:, 1:nch, g) = hub_n4jk(:, 1:nch, g) + part_n4(:, :, ig)
      hub_r4jk(:, 1:nch, g) = hub_r4jk(:, 1:nch, g) + part_r4(:, :, ig)
    end do
  end subroutine fold_4pcf_partials

end module query_4pcf_gpu_module
