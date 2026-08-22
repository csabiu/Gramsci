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
! Write routines and S4 symmetry table are reused from query_4pcf_module.
! ---------------------------------------------------------------------------
module query_4pcf_gpu_module
  use kdtree2_precision_module
  use iso_fortran_env, only: int8, int64
  use config_module
  use csr_module
  use query_4pcf_module, only: &
    write_4pcf_results, write_4pcf_results_noparity, &
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

    ! Runtime sizing
    integer :: max_nn, ngang, nwin, env_stat
    integer(int64) :: mnn2, lmat_bytes, freeb, win_edges
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

    ! ---- Runtime sizing -------------------------------------------------
    max_nn = max(csr_max_row_len(), 1)
    mnn2   = int(max_nn, int64)**2
    freeb  = gpu_free_mem_bytes()

    if (freeb < 0) then
      ngang = 1024        ! host fallback: no device limit, modest scratch
    else
      ! lmat scratch capped at ~1/3 of available device memory.
      ngang = int(min(4096_int64, max(512_int64, &
                      ((freeb - RESERVE_BYTES) / 3) / mnn2)))
    end if
    lmat_bytes = int(ngang, int64) * mnn2

    ! Three edge windows (hub + 2 search) resident in chunked mode.
    win_edges = csr_edge_window(3, lmat_bytes)
    nwin = int((csr_total_edges - 1) / win_edges) + 1
    ! Without an explicit override, prefer single-pass whenever the whole
    ! CSR fits alongside the scratch (= up to 3 windows' worth of edges).
    call get_environment_variable('GRAMSCI_GPU_WIN_EDGES', env, status=env_stat)
    if (env_stat /= 0 .and. nwin > 1 .and. nwin <= 3) nwin = 1

    if (cfg%rank == 0) &
      print '("4PCF GPU: max_nn=",i0,"  gangs=",i0,"  lmat scratch=",f6.2," GB",'// &
            '"  windows=",i0)', max_nn, ngang, real(lmat_bytes, kdkind)/1.0d9, nwin

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
      !$ACC& COPYIN(csr_ptr, csr_id, csr_dist, weights, buffer, bintable6) &
      !$ACC& CREATE(lmat_g) COPY(part_n4, part_r4)

      ! One gang per ig slot; hubs dealt round-robin so gang loads even out.
      ! Each gang exclusively owns lmat_g(:,:,ig) and part_*(:,ig).
      ! VECTOR_LENGTH(64): inner-loop trip counts shrink with k2/k3, so shorter
      ! vectors keep lanes busier than 128 while still filling the SMs.
      !$ACC PARALLEL LOOP GANG VECTOR_LENGTH(64) &
      !$ACC& PRIVATE(ig, i, k1, k2, nn_i, id1, id2, base_i, &
      !$ACC&         ind1, ind2, ind4, w12, rand12)
      do ig = 1, ngang
        do i = istart + ig - 1, iend, ngang
          if (buffer(i) == 1) cycle
          nn_i = int(csr_ptr(i + 1) - csr_ptr(i))
          if (nn_i <= 2) cycle

          base_i = csr_ptr(i) - 1

          ! ---- Phase 1: precompute local connectivity matrix ---------------
          ! C(nn_i, 2) binary searches.  For fixed k1 all lanes search the SAME
          ! id1 adjacency list (different targets) → warp-coherent, L1-cached.
          do k1 = 1, nn_i
            id1 = csr_id(base_i + k1)
            !$ACC LOOP VECTOR PRIVATE(k2, id2)
            do k2 = k1 + 1, nn_i
              id2 = csr_id(base_i + k2)
              lmat_g(k2, k1, ig) = &
                csr_find_dist_bin(csr_ptr, csr_id, csr_dist, id1, id2)
            end do
          end do

          ! ---- Phase 2: triple loop; k3 across lanes, lmat reads coalesced -
          do k1 = 1, nn_i
            ind1 = csr_dist(base_i + k1)
            id1  = csr_id  (base_i + k1)

            do k2 = k1 + 1, nn_i
              ind4 = lmat_g(k2, k1, ig)   ! O(1) lookup: dist bin of k1→k2
              if (ind4 == 0) cycle
              ind2 = csr_dist(base_i + k2)
              id2  = csr_id  (base_i + k2)
              w12  = weights(i) * weights(id1) * weights(id2)
              rand12 = (i > num_data_g .and. id1 > num_data_g .and. &
                        id2 > num_data_g)

              !$ACC LOOP VECTOR PRIVATE(k3, id3, config_idx, ind3, ind5, ind6, w4)
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
      !$ACC& COPYIN(csr_ptr, weights, buffer, bintable6) &
      !$ACC& CREATE(lmat_g) COPY(part_n4, part_r4)

      do cw1 = 1, nwin
        w1lo = split_id(cw1)
        w1hi = split_id(cw1 + 1) - 1
        off1 = split_e(cw1) - 1
        n1   = split_e(cw1 + 1) - split_e(cw1)
        stage1_id  (1:n1) = csr_id  (split_e(cw1):split_e(cw1 + 1) - 1)
        stage1_dist(1:n1) = csr_dist(split_e(cw1):split_e(cw1 + 1) - 1)

        !$ACC DATA COPYIN(stage1_id(1:n1), stage1_dist(1:n1))
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
          do hw = 1, nwin
            hublo = max(istart, split_id(hw))
            hubhi = min(iend,   split_id(hw + 1) - 1)
            if (hublo > hubhi) cycle
            elo_h = split_e(hw)
            ehi_h = split_e(hw + 1) - 1

            !$ACC DATA COPYIN(csr_id(elo_h:ehi_h), csr_dist(elo_h:ehi_h))
            !$ACC PARALLEL LOOP GANG VECTOR_LENGTH(64) &
            !$ACC& PRIVATE(ig, i, k1, k2, nn_i, id1, id2, base_i, &
            !$ACC&         ind1, ind2, ind4, w12, rand12)
            do ig = 1, ngang
              do i = hublo + ig - 1, hubhi, ngang
                if (buffer(i) == 1) cycle
                nn_i = int(csr_ptr(i + 1) - csr_ptr(i))
                if (nn_i <= 2) cycle

                base_i = csr_ptr(i) - 1

                ! ---- Phase 1a: lmat columns for neighbors in window cw1 ----
                do k1 = 1, nn_i
                  id1 = csr_id(base_i + k1)
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

                ! ---- Phase 2: pairs with id1 in cw1, id2 in cw2 ------------
                do k1 = 1, nn_i
                  id1 = csr_id(base_i + k1)
                  if (id1 < w1lo .or. id1 > w1hi) cycle
                  ind1 = csr_dist(base_i + k1)

                  do k2 = k1 + 1, nn_i
                    id2 = csr_id(base_i + k2)
                    if (id2 < w2lo .or. id2 > w2hi) cycle
                    ind4 = lmat_g(k2, k1, ig)
                    if (ind4 == 0) cycle
                    ind2 = csr_dist(base_i + k2)
                    w12  = weights(i) * weights(id1) * weights(id2)
                    rand12 = (i > num_data_g .and. id1 > num_data_g .and. &
                              id2 > num_data_g)

                    !$ACC LOOP VECTOR PRIVATE(k3, id3, config_idx, ind3, ind5, ind6, w4)
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

                    end do  ! k3 (VECTOR)
                  end do  ! k2
                end do  ! k1

              end do  ! i  (gang-sequential hubs)
            end do  ! ig  (GANG)
            !$ACC END PARALLEL LOOP
            !$ACC END DATA
          end do  ! hw
          !$ACC END DATA
        end do  ! cw2
        !$ACC END DATA
      end do  ! cw1
      !$ACC END DATA

      deallocate(stage1_id, stage1_dist, stage2_id, stage2_dist)
      deallocate(split_id, split_e)
    end if

    do config_idx = 1, n_configs_g
      N4(config_idx, 1) = N4(config_idx, 1) + sum(part_n4(config_idx, :))
      R4(config_idx, 1) = R4(config_idx, 1) + sum(part_r4(config_idx, :))
    end do

    deallocate(lmat_g, part_n4, part_r4)

    call write_4pcf_results_noparity()
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
    real(kdkind) :: w12, w4, vol
    real(kdkind) :: u1x, u1y, u1z, u2x, u2y, u2z, u3x, u3y, u3z, rn1, rn3
    logical :: rand12, exact_g
    integer :: num_data_g, n_configs_g

    ! Runtime sizing
    integer :: max_nn, ngang, nwin, env_stat
    integer(int64) :: mnn2, lmat_bytes, freeb, win_edges
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
    ! dummies so the COPYIN clauses are always valid
    if (.not. allocated(px)) allocate(px(1), py(1), pz(1))
    if (.not. allocated(csr_phi)) allocate(csr_phi(1))

    ! ---- Runtime sizing (see query_graph_4pcf_gpu) ----------------------
    max_nn = max(csr_max_row_len(), 1)
    mnn2   = int(max_nn, int64)**2
    freeb  = gpu_free_mem_bytes()

    if (freeb < 0) then
      ngang = 1024
    else
      ngang = int(min(4096_int64, max(512_int64, &
                      ((freeb - RESERVE_BYTES) / 3) / mnn2)))
    end if
    lmat_bytes = int(ngang, int64) * mnn2

    ! 4 window budgets: hub window costs 7 B/edge (id 4 + dist 1 + phi 2), the
    ! two staged search windows 5 B/edge — a 4×5 B budget covers 7+5+5 with
    ! slack.
    win_edges = csr_edge_window(4, lmat_bytes)
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
    allocate(part_n4(n_configs_g, 2, ngang))
    allocate(part_r4(n_configs_g, 2, ngang))
    part_n4 = 0.0d0
    part_r4 = 0.0d0

    if (nwin == 1) then
      ! =====================================================================
      ! Single-pass path
      ! =====================================================================
      !$ACC DATA &
      !$ACC& COPYIN(csr_ptr, csr_id, csr_dist, csr_phi, weights, buffer, &
      !$ACC&        bintable6, dir_x, dir_y, dir_z, px, py, pz) &
      !$ACC& CREATE(lmat_g, ug) COPY(part_n4, part_r4)

      !$ACC PARALLEL LOOP GANG VECTOR_LENGTH(64) &
      !$ACC& PRIVATE(ig, i, k1, k2, nn_i, id1, id2, base_i, &
      !$ACC&         ind1, ind2, ind4, p1, p2, w12, rand12, &
      !$ACC&         u1x, u1y, u1z, u2x, u2y, u2z, rn1)
      do ig = 1, ngang
        do i = istart + ig - 1, iend, ngang
          if (buffer(i) == 1) cycle
          nn_i = int(csr_ptr(i + 1) - csr_ptr(i))
          if (nn_i <= 2) cycle

          base_i = csr_ptr(i) - 1

          ! ---- Phase 1: connectivity matrix ------------------------------
          do k1 = 1, nn_i
            id1 = csr_id(base_i + k1)
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

          ! ---- Phase 2: triple loop with parity channel ------------------
          do k1 = 1, nn_i
            ind1 = csr_dist(base_i + k1)
            id1  = csr_id  (base_i + k1)
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
              !$ACC&   u3x, u3y, u3z, rn3)
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
                parity_flip = sign(1, raw_bin)

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

                !$ACC ATOMIC UPDATE
                part_n4(config_idx, 1, ig) = part_n4(config_idx, 1, ig) + w4
                !$ACC ATOMIC UPDATE
                part_n4(config_idx, 2, ig) = part_n4(config_idx, 2, ig) &
                                             + (parity_flip * sign_V) * w4

                if (rand12 .and. id3 > num_data_g) then
                  !$ACC ATOMIC UPDATE
                  part_r4(config_idx, 1, ig) = part_r4(config_idx, 1, ig) + w4
                  !$ACC ATOMIC UPDATE
                  part_r4(config_idx, 2, ig) = part_r4(config_idx, 2, ig) &
                                               + (parity_flip * sign_V) * w4
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
      !$ACC& COPYIN(csr_ptr, weights, buffer, bintable6, dir_x, dir_y, dir_z, &
      !$ACC&        px, py, pz) &
      !$ACC& CREATE(lmat_g, ug) COPY(part_n4, part_r4)

      do cw1 = 1, nwin
        w1lo = split_id(cw1)
        w1hi = split_id(cw1 + 1) - 1
        off1 = split_e(cw1) - 1
        n1   = split_e(cw1 + 1) - split_e(cw1)
        stage1_id  (1:n1) = csr_id  (split_e(cw1):split_e(cw1 + 1) - 1)
        stage1_dist(1:n1) = csr_dist(split_e(cw1):split_e(cw1 + 1) - 1)

        !$ACC DATA COPYIN(stage1_id(1:n1), stage1_dist(1:n1))
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
          do hw = 1, nwin
            hublo = max(istart, split_id(hw))
            hubhi = min(iend,   split_id(hw + 1) - 1)
            if (hublo > hubhi) cycle
            elo_h = split_e(hw)
            ehi_h = split_e(hw + 1) - 1
            ! csr_phi is a 1-element dummy under -exactparity (not stored)
            if (exact_g) then
              pelo = 1_int64
              pehi = 1_int64
            else
              pelo = elo_h
              pehi = ehi_h
            end if

            !$ACC DATA COPYIN(csr_id(elo_h:ehi_h), csr_dist(elo_h:ehi_h), &
            !$ACC&            csr_phi(pelo:pehi))
            !$ACC PARALLEL LOOP GANG VECTOR_LENGTH(64) &
            !$ACC& PRIVATE(ig, i, k1, k2, nn_i, id1, id2, base_i, &
            !$ACC&         ind1, ind2, ind4, p1, p2, w12, rand12, &
            !$ACC&         u1x, u1y, u1z, u2x, u2y, u2z, rn1)
            do ig = 1, ngang
              do i = hublo + ig - 1, hubhi, ngang
                if (buffer(i) == 1) cycle
                nn_i = int(csr_ptr(i + 1) - csr_ptr(i))
                if (nn_i <= 2) cycle

                base_i = csr_ptr(i) - 1

                ! ---- Phase 1a: columns for neighbors in window cw1 ---------
                do k1 = 1, nn_i
                  id1 = csr_id(base_i + k1)
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

                ! ---- Phase 2: pairs with id1 in cw1, id2 in cw2 ------------
                do k1 = 1, nn_i
                  id1 = csr_id(base_i + k1)
                  if (id1 < w1lo .or. id1 > w1hi) cycle
                  ind1 = csr_dist(base_i + k1)
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
                    !$ACC&   u3x, u3y, u3z, rn3)
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
                      parity_flip = sign(1, raw_bin)

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

                      !$ACC ATOMIC UPDATE
                      part_n4(config_idx, 1, ig) = part_n4(config_idx, 1, ig) + w4
                      !$ACC ATOMIC UPDATE
                      part_n4(config_idx, 2, ig) = part_n4(config_idx, 2, ig) &
                                                   + (parity_flip * sign_V) * w4

                      if (rand12 .and. id3 > num_data_g) then
                        !$ACC ATOMIC UPDATE
                        part_r4(config_idx, 1, ig) = part_r4(config_idx, 1, ig) + w4
                        !$ACC ATOMIC UPDATE
                        part_r4(config_idx, 2, ig) = part_r4(config_idx, 2, ig) &
                                                     + (parity_flip * sign_V) * w4
                      end if

                    end do  ! k3 (VECTOR)
                  end do  ! k2
                end do  ! k1

              end do  ! i
            end do  ! ig (GANG)
            !$ACC END PARALLEL LOOP
            !$ACC END DATA
          end do  ! hw
          !$ACC END DATA
        end do  ! cw2
        !$ACC END DATA
      end do  ! cw1
      !$ACC END DATA

      deallocate(stage1_id, stage1_dist, stage2_id, stage2_dist)
      deallocate(split_id, split_e)
    end if

    do config_idx = 1, n_configs_g
      N4(config_idx, 1) = N4(config_idx, 1) + sum(part_n4(config_idx, 1, :))
      N4(config_idx, 2) = N4(config_idx, 2) + sum(part_n4(config_idx, 2, :))
      R4(config_idx, 1) = R4(config_idx, 1) + sum(part_r4(config_idx, 1, :))
      R4(config_idx, 2) = R4(config_idx, 2) + sum(part_r4(config_idx, 2, :))
    end do

    deallocate(lmat_g, part_n4, part_r4)
    if (allocated(ug)) deallocate(ug)

    call write_4pcf_results()
  end subroutine query_graph_4pcf_parity_gpu

end module query_4pcf_gpu_module
