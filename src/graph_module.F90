module graph_module
  use kdtree2_module
  use kdtree2_precision_module
  use node_module
  use sorting_module
  use iso_fortran_env, only: int8, int16
  use config_module
  implicit none

  real(kdkind), parameter :: PI = 3.141592653589793238462643383279502884197d0

  ! Module-level KD-tree pointer (set by main before calling create_graph)
  type(kdtree2), pointer :: active_kd_tree => null()

contains

  subroutine set_kd_tree(tree)
    type(kdtree2), pointer, intent(in) :: tree
    active_kd_tree => tree
  end subroutine set_kd_tree

  subroutine create_graph(istart, iend)
    integer, intent(in) :: istart, iend

    integer :: i, j, k, nn1, nn2, nnode
    real(kdkind) :: vec_dist, theta
    real(kdkind) :: neighbor_pt(3), current_pt(3), mid_pt(3), work_vec(3)
    real(kdkind) :: dx, dy, dz, theta_sph, phi_sph
    integer :: i_theta, i_phi
    integer(int8) :: mu1_local
    logical :: store_phi
    type(kdtree2_result), allocatable :: resultsb(:)
    ! Periodic mode: which image-query center each kept neighbor came from
    integer(int8), allocatable :: rescen(:)
    integer :: c, ncenters, nf, sx, sy, sz, ntotal, iend_clamped, cap0
    real(kdkind) :: centers(3, 8), shift(3)
    integer :: bin_m, jbin
    real(kdkind) :: mu2, wpair

    ntotal = cfg%num_data + cfg%num_rand
    ! The drivers call create_graph(1, 999) for the timing estimate even when
    ! fewer points exist; clamp so we never index beyond the catalogue.
    iend_clamped = min(iend, ntotal)
    ! Per-thread search scratch, sized from sampled neighbor counts and
    ! grown on demand inside the loop -- NOT O(ntotal).  resultsb is
    ! OpenMP-private (each thread's copy is allocated to these bounds on
    ! entering the parallel region), so an ntotal-sized buffer costs
    ! 16 B * N per thread: ~100 GB of scratch at N = 1e8 on 64 threads.
    cap0 = initial_neighbor_cap(ntotal)
    allocate(resultsb(cap0))
    if (cfg%periodic) allocate(rescen(cap0))
    ! -exactparity takes the signed volume from positions, so the per-edge
    ! direction index is not needed at all (saves 2 B/edge).
    store_phi = cfg%four_pcf_parity .and. .not. cfg%exact_parity

    ! Pair Legendre-multipole sums for the anisotropic disconnected-4PCF
    ! subtraction.  Accumulated here because this is the only place the
    ! full-precision pair mu exists; each hub is visited exactly once across
    ! the drivers' probe (1..999) and main (1000..) sweeps, so no reset.
    if (.not. allocated(sum_pair_l2)) then
      allocate(sum_pair_l2(cfg%nbins), sum_pair_l4(cfg%nbins))
      sum_pair_l2 = 0.0d0
      sum_pair_l4 = 0.0d0
    end if

    !$OMP PARALLEL DO schedule(dynamic) &
    !$OMP& private(i, j, k, vec_dist, neighbor_pt, current_pt, mid_pt, work_vec, &
    !$OMP&         nn1, nn2, nnode, resultsb, rescen, theta, mu1_local, &
    !$OMP&         dx, dy, dz, theta_sph, phi_sph, i_theta, i_phi, &
    !$OMP&         c, ncenters, nf, sx, sy, sz, centers, shift, bin_m, jbin, mu2, wpair) &
    !$OMP& shared(active_kd_tree, output, points, weights, buffer, cfg, store_phi, &
    !$OMP&        ntotal, iend_clamped) &
    !$OMP& reduction(+: sum_pair_l2, sum_pair_l4)
    do i = istart, iend_clamped

      current_pt = points(:, i)

      if (cfg%periodic) then
        ! Periodic minimum-image search: query the (unmodified, non-periodic)
        ! tree around the hub and around its periodic images near the box
        ! faces.  With rmax < L/2 each neighbor is found exactly once, and
        ! its displacement from the query center is the minimum-image one.
        ncenters = 1
        centers(:, 1) = current_pt
        do sz = merge(-1, 0, current_pt(3) > cfg%boxsize(3) - cfg%rmax), &
                merge( 1, 0, current_pt(3) < cfg%rmax)
        do sy = merge(-1, 0, current_pt(2) > cfg%boxsize(2) - cfg%rmax), &
                merge( 1, 0, current_pt(2) < cfg%rmax)
        do sx = merge(-1, 0, current_pt(1) > cfg%boxsize(1) - cfg%rmax), &
                merge( 1, 0, current_pt(1) < cfg%rmax)
          if (sx == 0 .and. sy == 0 .and. sz == 0) cycle
          shift = [real(sx, kdkind) * cfg%boxsize(1), &
                   real(sy, kdkind) * cfg%boxsize(2), &
                   real(sz, kdkind) * cfg%boxsize(3)]
          ncenters = ncenters + 1
          centers(:, ncenters) = current_pt + shift
        end do
        end do
        end do

        ! Aggregate + filter in place: after this loop resultsb(1:nnode)
        ! holds the kept neighbors and rescen their query-center index.
        nnode = 0
        do c = 1, ncenters
          call kdtree2_r_nearest(tp=active_kd_tree, qv=centers(:, c), &
               r2=cfg%rmax*cfg%rmax, nfound=nf, nalloc=size(resultsb)-nnode, &
               results=resultsb(nnode+1:))
          do while (nf > size(resultsb) - nnode)
            ! Scratch too small for this center: nf is the true count even
            ! on overflow, so grow (preserving the neighbors kept from the
            ! earlier centers) and redo this center's query.
            call grow_scratch(resultsb, rescen, nnode + nf + 64, nnode)
            call kdtree2_r_nearest(tp=active_kd_tree, qv=centers(:, c), &
                 r2=cfg%rmax*cfg%rmax, nfound=nf, nalloc=size(resultsb)-nnode, &
                 results=resultsb(nnode+1:))
          end do
          k = nnode
          do j = nnode + 1, nnode + nf
            if (resultsb(j)%idx == i) cycle
            if (cfg%half_graph .and. resultsb(j)%idx <= i) cycle
            if (resultsb(j)%dis <= cfg%rmin * cfg%rmin) cycle
            if (resultsb(j)%dis >= cfg%rmax * cfg%rmax) cycle
            k = k + 1
            resultsb(k) = resultsb(j)
            rescen(k) = int(c, int8)
          end do
          nnode = k
        end do
      else

      call kdtree2_r_nearest_around_point(tp=active_kd_tree, idxin=i, correltime=-1, &
           r2=cfg%rmax*cfg%rmax, nfound=nn2, nalloc=size(resultsb), &
           results=resultsb)
      do while (nn2 > size(resultsb))
        ! Scratch too small for this hub: nn2 is the true count even on
        ! overflow, so one exact re-allocation and retry suffices (the
        ! sorted-results guarantee holds on the clean final call).
        deallocate(resultsb)
        allocate(resultsb(nn2 + 64))
        call kdtree2_r_nearest_around_point(tp=active_kd_tree, idxin=i, correltime=-1, &
             r2=cfg%rmax*cfg%rmax, nfound=nn2, nalloc=size(resultsb), &
             results=resultsb)
      end do

      if (cfg%rmin > 0.0) then
        nn1 = kdtree2_r_count_around_point(tp=active_kd_tree, idxin=i, correltime=-1, &
              r2=cfg%rmin*cfg%rmin)
      else
        nn1 = 1
      end if
      ! Pass 1: count valid neighbors with half-graph filter applied.
      ! Uses squared-distance comparisons to avoid redundant sqrt calls.
      nnode = 0
      do j = nn1 + 1, nn2
        if (cfg%half_graph .and. resultsb(j)%idx <= i) cycle
        if (resultsb(j)%dis <= cfg%rmin * cfg%rmin) cycle
        if (resultsb(j)%dis >= cfg%rmax * cfg%rmax) cycle
        nnode = nnode + 1
      end do

      end if

      if (store_phi) then
        call output(i)%init_with_phi(nnode)
      else
        call output(i)%init(nnode)
      end if

      if (cfg%periodic) then
        ! Fill from the pre-filtered aggregate.  RSD is rejected with -box,
        ! so mu is always 1; direction pixels use the min-image displacement.
        do k = 1, nnode
          vec_dist = sqrt(resultsb(k)%dis)
          output(i)%mu(k) = 1
          ! The filter above compared squared distances; sqrt rounding can
          ! still land exactly on rmin/rmax, so clamp the bin into range
          ! (unlike the non-periodic pass, we cannot cycle here — the
          ! adjacency arrays were already sized to nnode).
          if (cfg%logbins) then
            j = int(floor((log10(vec_dist) - cfg%log_rmin) * cfg%inv_delta_r)) + 1
          else
            j = int(floor((vec_dist - cfg%rmin) * cfg%inv_delta_r)) + 1
          end if
          output(i)%dist(k) = int(min(cfg%nbins, max(1, j)), int8)
          output(i)%id(k) = resultsb(k)%idx

          if (cfg%disc_rsd .and. buffer(i) /= 1) then
            ! Plane-parallel line of sight (z axis): mu of the min-image pair.
            dz = points(3, resultsb(k)%idx) - centers(3, rescen(k))
            mu2 = dz * dz / resultsb(k)%dis
            wpair = weights(i) * weights(resultsb(k)%idx)
            bin_m = int(output(i)%dist(k))
            sum_pair_l2(bin_m) = sum_pair_l2(bin_m) &
              + wpair * (1.5d0 * mu2 - 0.5d0)
            sum_pair_l4(bin_m) = sum_pair_l4(bin_m) &
              + wpair * (4.375d0 * mu2 * mu2 - 3.75d0 * mu2 + 0.375d0)
          end if

          if (store_phi) then
            neighbor_pt = points(:, resultsb(k)%idx)
            dx = neighbor_pt(1) - centers(1, rescen(k))
            dy = neighbor_pt(2) - centers(2, rescen(k))
            dz = neighbor_pt(3) - centers(3, rescen(k))

            theta_sph = acos(max(-1.0d0, min(1.0d0, dz / vec_dist)))  ! [0, pi], clamped
            phi_sph = atan2(dy, dx)                                     ! [-pi, pi]

            i_theta = min(cfg%n_theta_dir, max(1, int(floor(theta_sph * cfg%n_theta_dir / PI)) + 1))
            i_phi   = min(cfg%n_phi_dir,   max(1, int(floor((phi_sph + PI) * cfg%n_phi_dir / (2.0d0*PI))) + 1))

            ! int16, matching node%phi and the non-periodic pass: pixel
            ! indices reach n_theta*n_phi (256 on the default 8x32 grid),
            ! and an int8 conversion wraps anything above 127.
            output(i)%phi(k) = int((i_theta - 1) * cfg%n_phi_dir + i_phi, int16)
          end if
        end do
        k = nnode
      else

      ! Pass 2: fill with EXACTLY the same squared-distance filter as pass 1
      ! (a sqrt-based filter here could disagree with pass 1 at values whose
      ! square rounds onto rmin^2/rmax^2 and overflow the arrays sized above).
      k = 0
      do j = nn1 + 1, nn2
        if (cfg%half_graph .and. resultsb(j)%idx <= i) cycle
        if (resultsb(j)%dis <= cfg%rmin * cfg%rmin) cycle
        if (resultsb(j)%dis >= cfg%rmax * cfg%rmax) cycle
        neighbor_pt = points(:, resultsb(j)%idx)

        vec_dist = sqrt(resultsb(j)%dis)
        k = k + 1

        if (cfg%RSD) then
          mid_pt = 0.5d0 * (current_pt + neighbor_pt)
          work_vec = current_pt - neighbor_pt
          theta = dot_product(mid_pt, work_vec) / &
                  (sqrt(dot_product(mid_pt, mid_pt)) * sqrt(dot_product(work_vec, work_vec)))

          mu1_local = int(floor((theta + 1.0d0) * cfg%mu_scale) + 1, int8)
          if (mu1_local < 1) mu1_local = 1
          if (mu1_local > cfg%nmu) mu1_local = int(cfg%nmu, int8)
          output(i)%mu(k) = mu1_local
        else
          output(i)%mu(k) = 1
        end if

        ! Clamp the bin like the periodic pass: the filter compared squared
        ! distances, so sqrt/reciprocal rounding can still push the index to
        ! 0 or nbins+1 at the bin-range edges.
        if (cfg%logbins) then
          jbin = int(floor((log10(vec_dist) - cfg%log_rmin) * cfg%inv_delta_r)) + 1
        else
          jbin = int(floor((vec_dist - cfg%rmin) * cfg%inv_delta_r)) + 1
        end if
        output(i)%dist(k) = int(min(cfg%nbins, max(1, jbin)), int8)
        output(i)%id(k) = resultsb(j)%idx

        ! Survey mode: disc_rsd implies cfg%RSD, so theta (the midpoint-LOS
        ! pair mu) was computed above.
        if (cfg%disc_rsd .and. buffer(i) /= 1) then
          mu2 = theta * theta
          wpair = weights(i) * weights(resultsb(j)%idx)
          bin_m = int(output(i)%dist(k))
          sum_pair_l2(bin_m) = sum_pair_l2(bin_m) &
            + wpair * (1.5d0 * mu2 - 0.5d0)
          sum_pair_l4(bin_m) = sum_pair_l4(bin_m) &
            + wpair * (4.375d0 * mu2 * mu2 - 3.75d0 * mu2 + 0.375d0)
        end if

        ! Compute and store direction pixel for 4PCF parity
        if (store_phi) then
          dx = neighbor_pt(1) - current_pt(1)
          dy = neighbor_pt(2) - current_pt(2)
          dz = neighbor_pt(3) - current_pt(3)

          ! Spherical coordinates in fixed Cartesian frame
          theta_sph = acos(max(-1.0d0, min(1.0d0, dz / vec_dist)))  ! [0, pi], clamped
          phi_sph = atan2(dy, dx)                                     ! [-pi, pi]

          ! Bin indices (1-based, clamped)
          i_theta = min(cfg%n_theta_dir, max(1, int(floor(theta_sph * cfg%n_theta_dir / PI)) + 1))
          i_phi   = min(cfg%n_phi_dir,   max(1, int(floor((phi_sph + PI) * cfg%n_phi_dir / (2.0d0*PI))) + 1))

          ! Combined pixel index (1 to n_dir_pixels)
          output(i)%phi(k) = int((i_theta - 1) * cfg%n_phi_dir + i_phi, int16)
        end if
      end do

      end if

      output(i)%nn = k
      if (k > 0) then
        if (store_phi) then
          call sort2_with_phi(output(i)%id, output(i)%dist, output(i)%mu, output(i)%phi, k)
        else
          call sort2(output(i)%id, output(i)%dist, output(i)%mu, k)
        end if
      end if

    end do
    !$OMP END PARALLEL DO

  end subroutine create_graph

  ! Sample hub neighbor counts to choose the initial per-thread search
  ! scratch size.  Twice the sampled maximum plus headroom absorbs density
  ! fluctuations and (in periodic mode) the image-center aggregation near
  ! the box faces; any hub that still exceeds it triggers the grow-and-
  ! retry path in create_graph, so the choice affects only performance.
  function initial_neighbor_cap(ntotal) result(cap)
    integer, intent(in) :: ntotal
    integer :: cap
    integer :: s, j, nc, nsample
    real(kdkind) :: rv

    nsample = min(200, ntotal)
    cap = 0
    do s = 1, nsample
      call random_number(rv)
      j = min(ntotal, int(rv * real(ntotal, kdkind)) + 1)
      nc = kdtree2_r_count_around_point(tp=active_kd_tree, idxin=j, &
           correltime=-1, r2=cfg%rmax*cfg%rmax)
      cap = max(cap, nc)
    end do
    cap = min(ntotal, 2 * cap + 64)
  end function initial_neighbor_cap

  ! Grow the per-thread search scratch to newcap entries, preserving the
  ! first nkeep (the filtered neighbors kept from earlier image centers).
  ! rescen is grown in lockstep so the two arrays never disagree in size;
  ! it is unallocated in non-periodic mode and then left untouched.
  subroutine grow_scratch(res, cen, newcap, nkeep)
    type(kdtree2_result), allocatable, intent(inout) :: res(:)
    integer(int8), allocatable, intent(inout) :: cen(:)
    integer, intent(in) :: newcap, nkeep
    type(kdtree2_result), allocatable :: rtmp(:)
    integer(int8), allocatable :: ctmp(:)

    allocate(rtmp(newcap))
    rtmp(1:nkeep) = res(1:nkeep)
    call move_alloc(rtmp, res)
    if (allocated(cen)) then
      allocate(ctmp(newcap))
      ctmp(1:nkeep) = cen(1:nkeep)
      call move_alloc(ctmp, cen)
    end if
  end subroutine grow_scratch

end module graph_module
