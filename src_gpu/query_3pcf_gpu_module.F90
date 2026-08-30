! ---------------------------------------------------------------------------
! GPU-offloaded 3PCF query using OpenACC.
!
! Strategy: binary-search (bsearch) over the CSR graph rather than the
! merge-walk used on CPU.  The merge-walk's sequential pointer state across a
! single triangle is hostile to warp execution; the O(log m) binary search
! overhead is negligible under GPU parallelism.
!
! Parallel mapping: one GANG (CUDA block) per hub node i; the inner k2 loop
! is spread across the gang's VECTOR lanes.  For a fixed k1 every lane
! binary-searches the SAME id1 adjacency list (for different targets), so the
! warp's memory traffic is confined to one small list segment → L1/L2-coherent.
! (Letting the compiler choose maps hubs onto gang AND vector, which puts a
! different hub on every lane: fully divergent accesses and warp-level load
! imbalance — measured 6× slower than CPU×64 at rmax=80.)
!
! Accumulation: vector lanes within a gang race on a shared accumulator, and
! nvfortran array reductions across gang+vector levels re-reduce at the end of
! every inner loop instance.  Instead, atomics into slot-strided partial
! arrays part_*(bin, slot) with slot = mod(i-1, NSLOT)+1: concurrently
! resident gangs (~1000) mostly own distinct slots, so atomic contention is
! negligible.  Slots are summed on the host afterwards.
!
! Chunked mode (graphs larger than device memory): the kernel needs two kinds
! of CSR row access — the hub's own row (sequential) and the row of any
! neighbor id1 (bsearch target).  The edge arrays are split into NWIN
! row-aligned windows sized at runtime from free device memory (override:
! GRAMSCI_GPU_WIN_EDGES) and the query runs as an
! NWIN × NWIN tile loop over (search-window, hub-window) pairs.  A triangle
! at hub i with pair (id1, id2) is processed exactly once: in the tile whose
! hub window holds i's row and whose search window holds id1's row.  Tiles
! skip a whole k1 iteration when id1 is outside the current search window
! (one integer compare), so duplicated work is negligible.  The search window
! is staged into separate host arrays because OpenACC cannot map two disjoint
! sections of the same array simultaneously.  Device footprint: 2 windows.
!
! Limitation: RSD mode (nmu > 1) falls back to the CPU implementation because
! find_normal's GPU equivalent is not yet implemented.
!
! Write routines are unchanged — they live in query_3pcf_module and run on CPU
! after the kernel completes and results are copied back to host.
! ---------------------------------------------------------------------------
module query_3pcf_gpu_module
  use kdtree2_precision_module
  use iso_fortran_env, only: int8, int64
  use config_module
  use csr_module
  use query_3pcf_module, only: write_3pcf_results, write_equilateral_results, &
                               write_3pcf_jackknife
  implicit none

  ! Number of independent accumulator copies; must comfortably exceed the
  ! count of concurrently resident gangs (16 blocks/SM × SM count; 4096
  ! covers >250-SM parts, and the host-side cost of extra slots is trivial).
  integer, parameter :: NSLOT_3PCF = 4096
  ! The jackknife partials are indexed (bin, slot, region), so they cannot use
  ! NSLOT_3PCF: at 84 bins x 400 regions that would be 1.1 TB. The region index
  ! already spreads atomics ~njk ways, so far fewer slots suffice.
  integer, parameter :: NSLOT_JK = 64
  ! Refuse rather than silently exhaust device memory.
  real(kdkind), parameter :: JK_MAX_GB = 4.0d0

contains

  subroutine query_graph_3pcf_gpu(istart, iend)
    integer, intent(in) :: istart, iend

    integer :: i, k1, k2, nn_i, id1, id2, bin, slot
    integer(int64) :: base_i
    integer(int8) :: ind1, ind2, ind3
    real(kdkind) :: wi_w1, w3
    integer :: num_data_g, config_bins_g

    ! Chunking bookkeeping
    integer :: nwin, cw, hw
    integer :: idlo2, idhi2, hublo, hubhi
    integer(int64) :: off2, elo_h, ehi_h, nstage, maxw, win_edges
    integer, allocatable :: split_id(:)
    integer(int64), allocatable :: split_e(:)
    integer, allocatable :: stage_id(:)
    integer(int8), allocatable :: stage_dist(:)

    ! Slot-strided partial accumulators (config_bins, NSLOT_3PCF); summed
    ! over slots on the host after the kernel.
    real(kdkind), allocatable :: part_nnn(:,:), part_rrr(:,:)
    ! Jackknife partials (bin, slot, region) and per-triplet region labels.
    real(kdkind), allocatable :: part_nnn_jk(:,:,:), part_rrr_jk(:,:,:)
    integer :: njk_g, jslot, jr1, jr2, jr3, m
    real(kdkind) :: jk_gb

    if (cfg%RSD) then
      print *, 'ERROR: query_graph_3pcf_gpu called with RSD mode; caller must route to CPU'
      stop 1
    end if

    if (cfg%rank == 0) print *, 'Performing 3PCF (all configurations, GPU bsearch)'
    if (cfg%rank == 0) print *, 'begin querying the graph'

    num_data_g    = cfg%num_data
    config_bins_g = cfg%config_bins
    njk_g         = cfg%njk

    if (njk_g > 0) then
      jk_gb = 2.0d0 * real(config_bins_g, kdkind) * real(NSLOT_JK, kdkind) * &
              real(njk_g, kdkind) * 8.0d0 / 1.0d9
      if (cfg%rank == 0) print *, 'jackknife partials: ', jk_gb, ' GB on device'
      if (jk_gb > JK_MAX_GB) then
        print *, 'ERROR: jackknife partial arrays need ', jk_gb, &
                 ' GB (limit ', JK_MAX_GB, '). Reduce -njk or the bin count, '// &
                 'or lower NSLOT_JK.'
        stop
      end if
    end if
    ! Allocated unconditionally: the ACC data clauses name them either way.
    allocate(part_nnn_jk(config_bins_g, NSLOT_JK, max(njk_g, 1)))
    allocate(part_rrr_jk(config_bins_g, NSLOT_JK, max(njk_g, 1)))
    part_nnn_jk = 0.0d0
    part_rrr_jk = 0.0d0

    allocate(part_nnn(config_bins_g, NSLOT_3PCF))
    allocate(part_rrr(config_bins_g, NSLOT_3PCF))
    part_nnn = 0.0d0
    part_rrr = 0.0d0

    ! Two edge windows (hub + search) resident in chunked mode.
    win_edges = csr_edge_window(2, 0_int64)
    nwin = int((csr_total_edges - 1) / win_edges) + 1

    if (nwin == 1) then
      ! =====================================================================
      ! Single-pass path: whole CSR fits on the device.
      ! =====================================================================
      !$ACC DATA &
      !$ACC& COPYIN(csr_ptr, csr_id, csr_dist, weights, bintable, region) &
      !$ACC& COPY(part_nnn, part_rrr, part_nnn_jk, part_rrr_jk)

      !$ACC PARALLEL LOOP GANG VECTOR_LENGTH(128) &
      !$ACC& PRIVATE(i, k1, nn_i, id1, base_i, slot, ind1, wi_w1, jslot, jr1)
      do i = istart, iend

        nn_i = int(csr_ptr(i + 1) - csr_ptr(i))
        if (nn_i < 2) cycle

        base_i = csr_ptr(i) - 1
        slot   = mod(i - 1, NSLOT_3PCF) + 1
        jslot  = mod(i - 1, NSLOT_JK) + 1
        if (njk_g > 0) jr1 = region(i)

        do k1 = 1, nn_i
          ind1  = csr_dist(base_i + k1)
          id1   = csr_id  (base_i + k1)
          wi_w1 = weights(i) * weights(id1)

          ! All lanes search id1's adjacency list — same segment, coherent.
          !$ACC LOOP VECTOR PRIVATE(k2, id2, bin, ind2, ind3, w3, jr2, jr3)
          do k2 = k1 + 1, nn_i
            ind2 = csr_dist(base_i + k2)
            id2  = csr_id  (base_i + k2)

            ind3 = csr_find_dist_bin(csr_ptr, csr_id, csr_dist, id1, id2)
            if (ind3 == 0) cycle

            bin = bintable(int(ind1), int(ind2), int(ind3), 1)
            w3  = wi_w1 * weights(id2)

            !$ACC ATOMIC UPDATE
            part_nnn(bin, slot) = part_nnn(bin, slot) + w3

            if (njk_g > 0) then
              jr2 = region(id1)
              jr3 = region(id2)
              if (jr1 > 0) then
                !$ACC ATOMIC UPDATE
                part_nnn_jk(bin, jslot, jr1) = part_nnn_jk(bin, jslot, jr1) + w3
              end if
              if (jr2 > 0 .and. jr2 /= jr1) then
                !$ACC ATOMIC UPDATE
                part_nnn_jk(bin, jslot, jr2) = part_nnn_jk(bin, jslot, jr2) + w3
              end if
              if (jr3 > 0 .and. jr3 /= jr1 .and. jr3 /= jr2) then
                !$ACC ATOMIC UPDATE
                part_nnn_jk(bin, jslot, jr3) = part_nnn_jk(bin, jslot, jr3) + w3
              end if
            end if

            if (i > num_data_g .and. id1 > num_data_g .and. id2 > num_data_g) then
              !$ACC ATOMIC UPDATE
              part_rrr(bin, slot) = part_rrr(bin, slot) - w3
              if (njk_g > 0) then
                if (jr1 > 0) then
                  !$ACC ATOMIC UPDATE
                  part_rrr_jk(bin, jslot, jr1) = part_rrr_jk(bin, jslot, jr1) - w3
                end if
                if (jr2 > 0 .and. jr2 /= jr1) then
                  !$ACC ATOMIC UPDATE
                  part_rrr_jk(bin, jslot, jr2) = part_rrr_jk(bin, jslot, jr2) - w3
                end if
                if (jr3 > 0 .and. jr3 /= jr1 .and. jr3 /= jr2) then
                  !$ACC ATOMIC UPDATE
                  part_rrr_jk(bin, jslot, jr3) = part_rrr_jk(bin, jslot, jr3) - w3
                end if
              end if
            end if
          end do
        end do
      end do
      !$ACC END PARALLEL LOOP
      !$ACC END DATA

    else
      ! =====================================================================
      ! Chunked path: NWIN × NWIN tiles, 2 edge windows resident at a time.
      ! =====================================================================
      call csr_make_splits(win_edges, nwin, split_id, split_e)

      maxw = 0
      do cw = 1, nwin
        maxw = max(maxw, split_e(cw + 1) - split_e(cw))
      end do
      if (cfg%rank == 0) &
        print '("3PCF GPU chunked: ",i0,"x",i0," tiles, window <= ",f6.2," GB")', &
              nwin, nwin, real(maxw, kdkind) * 5.0d0 / 1.0d9

      allocate(stage_id(maxw), stage_dist(maxw))

      !$ACC DATA &
      !$ACC& COPYIN(csr_ptr, weights, bintable, region) &
      !$ACC& COPY(part_nnn, part_rrr, part_nnn_jk, part_rrr_jk)

      do cw = 1, nwin
        ! ---- Stage the search window (rows whose lists get bsearched) ----
        idlo2 = split_id(cw)
        idhi2 = split_id(cw + 1) - 1
        off2  = split_e(cw) - 1
        nstage = split_e(cw + 1) - split_e(cw)
        stage_id  (1:nstage) = csr_id  (split_e(cw):split_e(cw + 1) - 1)
        stage_dist(1:nstage) = csr_dist(split_e(cw):split_e(cw + 1) - 1)

        !$ACC DATA COPYIN(stage_id(1:nstage), stage_dist(1:nstage))
        do hw = 1, nwin
          ! ---- Map the hub window (rows iterated as hubs) ----
          hublo = max(istart, split_id(hw))
          hubhi = min(iend,   split_id(hw + 1) - 1)
          if (hublo > hubhi) cycle
          elo_h = split_e(hw)
          ehi_h = split_e(hw + 1) - 1

          !$ACC DATA COPYIN(csr_id(elo_h:ehi_h), csr_dist(elo_h:ehi_h))
          !$ACC PARALLEL LOOP GANG VECTOR_LENGTH(128) &
          !$ACC& PRIVATE(i, k1, nn_i, id1, base_i, slot, ind1, wi_w1, jslot, jr1)
          do i = hublo, hubhi

            nn_i = int(csr_ptr(i + 1) - csr_ptr(i))
            if (nn_i < 2) cycle

            base_i = csr_ptr(i) - 1
            slot   = mod(i - 1, NSLOT_3PCF) + 1
            jslot  = mod(i - 1, NSLOT_JK) + 1
            if (njk_g > 0) jr1 = region(i)

            do k1 = 1, nn_i
              id1 = csr_id(base_i + k1)
              ! Pairs with this id1 belong to the tile holding id1's row.
              if (id1 < idlo2 .or. id1 > idhi2) cycle

              ind1  = csr_dist(base_i + k1)
              wi_w1 = weights(i) * weights(id1)

              !$ACC LOOP VECTOR PRIVATE(k2, id2, bin, ind2, ind3, w3, jr2, jr3)
              do k2 = k1 + 1, nn_i
                ind2 = csr_dist(base_i + k2)
                id2  = csr_id  (base_i + k2)

                ind3 = csr_find_dist_bin_win(csr_ptr, stage_id, stage_dist, &
                                             off2, id1, id2)
                if (ind3 == 0) cycle

                bin = bintable(int(ind1), int(ind2), int(ind3), 1)
                w3  = wi_w1 * weights(id2)

                !$ACC ATOMIC UPDATE
                part_nnn(bin, slot) = part_nnn(bin, slot) + w3

                if (njk_g > 0) then
                  jr2 = region(id1)
                  jr3 = region(id2)
                  if (jr1 > 0) then
                    !$ACC ATOMIC UPDATE
                    part_nnn_jk(bin, jslot, jr1) = part_nnn_jk(bin, jslot, jr1) + w3
                  end if
                  if (jr2 > 0 .and. jr2 /= jr1) then
                    !$ACC ATOMIC UPDATE
                    part_nnn_jk(bin, jslot, jr2) = part_nnn_jk(bin, jslot, jr2) + w3
                  end if
                  if (jr3 > 0 .and. jr3 /= jr1 .and. jr3 /= jr2) then
                    !$ACC ATOMIC UPDATE
                    part_nnn_jk(bin, jslot, jr3) = part_nnn_jk(bin, jslot, jr3) + w3
                  end if
                end if

                if (i > num_data_g .and. id1 > num_data_g .and. &
                    id2 > num_data_g) then
                  !$ACC ATOMIC UPDATE
                  part_rrr(bin, slot) = part_rrr(bin, slot) - w3
                  if (njk_g > 0) then
                    if (jr1 > 0) then
                      !$ACC ATOMIC UPDATE
                      part_rrr_jk(bin, jslot, jr1) = part_rrr_jk(bin, jslot, jr1) - w3
                    end if
                    if (jr2 > 0 .and. jr2 /= jr1) then
                      !$ACC ATOMIC UPDATE
                      part_rrr_jk(bin, jslot, jr2) = part_rrr_jk(bin, jslot, jr2) - w3
                    end if
                    if (jr3 > 0 .and. jr3 /= jr1 .and. jr3 /= jr2) then
                      !$ACC ATOMIC UPDATE
                      part_rrr_jk(bin, jslot, jr3) = part_rrr_jk(bin, jslot, jr3) - w3
                    end if
                  end if
                end if
              end do
            end do
          end do
          !$ACC END PARALLEL LOOP
          !$ACC END DATA
        end do  ! hw
        !$ACC END DATA
      end do  ! cw
      !$ACC END DATA

      deallocate(stage_id, stage_dist, split_id, split_e)
    end if

    ! Fold slot partials into N2/N3 (slot (:,1,3) = NNN / RRR counts).
    do bin = 1, config_bins_g
      N2(bin, 1, 3) = N2(bin, 1, 3) + sum(part_nnn(bin, :))
      N3(bin, 1, 3) = N3(bin, 1, 3) + sum(part_rrr(bin, :))
    end do
    if (njk_g > 0) then
      do m = 1, njk_g
        do bin = 1, config_bins_g
          N2jk(bin, 1, 3, m) = N2jk(bin, 1, 3, m) + sum(part_nnn_jk(bin, :, m))
          N3jk(bin, 1, 3, m) = N3jk(bin, 1, 3, m) + sum(part_rrr_jk(bin, :, m))
        end do
      end do
    end if

    deallocate(part_nnn, part_rrr, part_nnn_jk, part_rrr_jk)

    call write_3pcf_results()
    call write_3pcf_jackknife()
  end subroutine query_graph_3pcf_gpu


  ! ---------------------------------------------------------------------------
  ! Equilateral 3PCF on the GPU: same gang/vector mapping and chunked tiling
  ! as query_graph_3pcf_gpu, but pairs are filtered to ind1 == ind2 BEFORE the
  ! binary search (only ~1/nbins of pairs get searched) and the triangle is
  ! kept only when the third side falls in the same bin.  Counts accumulate
  ! per radial bin (no bintable lookup); output via write_equilateral_results.
  ! RSD mode must be routed to the CPU implementation by the caller.
  ! ---------------------------------------------------------------------------
  subroutine query_graph_equilateral_gpu(istart, iend)
    integer, intent(in) :: istart, iend
    ! -njk is implemented for the full 3PCF kernel only; refuse rather than
    ! return zero-scatter realisations.

    integer :: i, k1, k2, nn_i, id1, id2, bin, slot
    integer(int64) :: base_i
    integer(int8) :: ind1, ind2, ind3
    real(kdkind) :: wi_w1, w3
    integer :: num_data_g, nbins_g

    integer :: nwin, cw, hw
    integer :: idlo2, idhi2, hublo, hubhi
    integer(int64) :: off2, elo_h, ehi_h, nstage, maxw, win_edges
    integer, allocatable :: split_id(:)
    integer(int64), allocatable :: split_e(:)
    integer, allocatable :: stage_id(:)
    integer(int8), allocatable :: stage_dist(:)

    real(kdkind), allocatable :: part_nnn(:,:), part_rrr(:,:)

    if (cfg%RSD) then
      print *, 'ERROR: query_graph_equilateral_gpu called with RSD mode; caller must route to CPU'
      stop 1
    end if
    if (cfg%njk > 0) then
      print *, 'ERROR: -njk is implemented for -3pcf, not -equi.'
      stop
    end if

    if (cfg%rank == 0) print *, 'Performing equilateral 3PCF (GPU bsearch)'
    if (cfg%rank == 0) print *, 'begin querying the graph'

    num_data_g = cfg%num_data
    nbins_g    = cfg%nbins

    allocate(part_nnn(nbins_g, NSLOT_3PCF))
    allocate(part_rrr(nbins_g, NSLOT_3PCF))
    part_nnn = 0.0d0
    part_rrr = 0.0d0

    ! Two edge windows (hub + search) resident in chunked mode.
    win_edges = csr_edge_window(2, 0_int64)
    nwin = int((csr_total_edges - 1) / win_edges) + 1

    if (nwin == 1) then
      ! =====================================================================
      ! Single-pass path: whole CSR fits on the device.
      ! =====================================================================
      !$ACC DATA &
      !$ACC& COPYIN(csr_ptr, csr_id, csr_dist, weights) &
      !$ACC& COPY(part_nnn, part_rrr)

      !$ACC PARALLEL LOOP GANG VECTOR_LENGTH(128) &
      !$ACC& PRIVATE(i, k1, nn_i, id1, base_i, slot, ind1, wi_w1)
      do i = istart, iend

        nn_i = int(csr_ptr(i + 1) - csr_ptr(i))
        if (nn_i < 2) cycle

        base_i = csr_ptr(i) - 1
        slot   = mod(i - 1, NSLOT_3PCF) + 1

        do k1 = 1, nn_i
          ind1  = csr_dist(base_i + k1)
          id1   = csr_id  (base_i + k1)
          wi_w1 = weights(i) * weights(id1)

          !$ACC LOOP VECTOR PRIVATE(k2, id2, bin, ind2, ind3, w3)
          do k2 = k1 + 1, nn_i
            ind2 = csr_dist(base_i + k2)
            if (ind2 /= ind1) cycle      ! equilateral filter before bsearch

            id2  = csr_id(base_i + k2)
            ind3 = csr_find_dist_bin(csr_ptr, csr_id, csr_dist, id1, id2)
            if (ind3 /= ind1) cycle

            bin = int(ind1)
            w3  = wi_w1 * weights(id2)

            !$ACC ATOMIC UPDATE
            part_nnn(bin, slot) = part_nnn(bin, slot) + w3

            if (i > num_data_g .and. id1 > num_data_g .and. id2 > num_data_g) then
              !$ACC ATOMIC UPDATE
              part_rrr(bin, slot) = part_rrr(bin, slot) - w3
            end if
          end do
        end do
      end do
      !$ACC END PARALLEL LOOP
      !$ACC END DATA

    else
      ! =====================================================================
      ! Chunked path: NWIN × NWIN tiles, 2 edge windows resident at a time.
      ! =====================================================================
      call csr_make_splits(win_edges, nwin, split_id, split_e)

      maxw = 0
      do cw = 1, nwin
        maxw = max(maxw, split_e(cw + 1) - split_e(cw))
      end do
      if (cfg%rank == 0) &
        print '("equi 3PCF GPU chunked: ",i0,"x",i0," tiles, window <= ",f6.2," GB")', &
              nwin, nwin, real(maxw, kdkind) * 5.0d0 / 1.0d9

      allocate(stage_id(maxw), stage_dist(maxw))

      !$ACC DATA &
      !$ACC& COPYIN(csr_ptr, weights) &
      !$ACC& COPY(part_nnn, part_rrr)

      do cw = 1, nwin
        idlo2 = split_id(cw)
        idhi2 = split_id(cw + 1) - 1
        off2  = split_e(cw) - 1
        nstage = split_e(cw + 1) - split_e(cw)
        stage_id  (1:nstage) = csr_id  (split_e(cw):split_e(cw + 1) - 1)
        stage_dist(1:nstage) = csr_dist(split_e(cw):split_e(cw + 1) - 1)

        !$ACC DATA COPYIN(stage_id(1:nstage), stage_dist(1:nstage))
        do hw = 1, nwin
          hublo = max(istart, split_id(hw))
          hubhi = min(iend,   split_id(hw + 1) - 1)
          if (hublo > hubhi) cycle
          elo_h = split_e(hw)
          ehi_h = split_e(hw + 1) - 1

          !$ACC DATA COPYIN(csr_id(elo_h:ehi_h), csr_dist(elo_h:ehi_h))
          !$ACC PARALLEL LOOP GANG VECTOR_LENGTH(128) &
          !$ACC& PRIVATE(i, k1, nn_i, id1, base_i, slot, ind1, wi_w1)
          do i = hublo, hubhi

            nn_i = int(csr_ptr(i + 1) - csr_ptr(i))
            if (nn_i < 2) cycle

            base_i = csr_ptr(i) - 1
            slot   = mod(i - 1, NSLOT_3PCF) + 1

            do k1 = 1, nn_i
              id1 = csr_id(base_i + k1)
              if (id1 < idlo2 .or. id1 > idhi2) cycle

              ind1  = csr_dist(base_i + k1)
              wi_w1 = weights(i) * weights(id1)

              !$ACC LOOP VECTOR PRIVATE(k2, id2, bin, ind2, ind3, w3)
              do k2 = k1 + 1, nn_i
                ind2 = csr_dist(base_i + k2)
                if (ind2 /= ind1) cycle    ! equilateral filter before bsearch

                id2  = csr_id(base_i + k2)
                ind3 = csr_find_dist_bin_win(csr_ptr, stage_id, stage_dist, &
                                             off2, id1, id2)
                if (ind3 /= ind1) cycle

                bin = int(ind1)
                w3  = wi_w1 * weights(id2)

                !$ACC ATOMIC UPDATE
                part_nnn(bin, slot) = part_nnn(bin, slot) + w3

                if (i > num_data_g .and. id1 > num_data_g .and. &
                    id2 > num_data_g) then
                  !$ACC ATOMIC UPDATE
                  part_rrr(bin, slot) = part_rrr(bin, slot) - w3
                end if
              end do
            end do
          end do
          !$ACC END PARALLEL LOOP
          !$ACC END DATA
        end do  ! hw
        !$ACC END DATA
      end do  ! cw
      !$ACC END DATA

      deallocate(stage_id, stage_dist, split_id, split_e)
    end if

    ! Fold slot partials into N2/N3, indexed directly by radial bin.
    do bin = 1, nbins_g
      N2(bin, 1, 3) = N2(bin, 1, 3) + sum(part_nnn(bin, :))
      N3(bin, 1, 3) = N3(bin, 1, 3) + sum(part_rrr(bin, :))
    end do

    deallocate(part_nnn, part_rrr)

    call write_equilateral_results()
  end subroutine query_graph_equilateral_gpu

end module query_3pcf_gpu_module
