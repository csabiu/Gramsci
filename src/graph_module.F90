module graph_module
  use kdtree2_module
  use kdtree2_precision_module
  use node_module
  use sorting_module
  use iso_fortran_env, only: int8
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

    allocate(resultsb(cfg%num_data + cfg%num_rand))
    store_phi = cfg%four_pcf_parity

    !$OMP PARALLEL DO schedule(dynamic) &
    !$OMP& private(i, j, k, vec_dist, neighbor_pt, current_pt, mid_pt, work_vec, &
    !$OMP&         nn1, nn2, nnode, resultsb, theta, mu1_local, &
    !$OMP&         dx, dy, dz, theta_sph, phi_sph, i_theta, i_phi) &
    !$OMP& shared(active_kd_tree, output, points, cfg, store_phi)
    do i = istart, iend

      call kdtree2_r_nearest_around_point(tp=active_kd_tree, idxin=i, correltime=-1, &
           r2=cfg%rmax*cfg%rmax, nfound=nn2, nalloc=(cfg%num_data+cfg%num_rand), &
           results=resultsb)

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

      if (store_phi) then
        call output(i)%init_with_phi(nnode)
      else
        call output(i)%init(nnode)
      end if
      current_pt = points(:, i)

      ! Pass 2: fill with the same filter.
      k = 0
      do j = nn1 + 1, nn2
        if (cfg%half_graph .and. resultsb(j)%idx <= i) cycle
        neighbor_pt = points(:, resultsb(j)%idx)

        vec_dist = sqrt(resultsb(j)%dis)

        if (vec_dist <= cfg%rmin) cycle
        if (vec_dist >= cfg%rmax) cycle
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

        if (cfg%logbins) then
          output(i)%dist(k) = int(floor((log10(vec_dist) - cfg%log_rmin) * cfg%inv_delta_r) + 1, int8)
        else
          output(i)%dist(k) = int(floor((vec_dist - cfg%rmin) * cfg%inv_delta_r) + 1, int8)
        end if
        output(i)%id(k) = resultsb(j)%idx

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
          output(i)%phi(k) = int((i_theta - 1) * cfg%n_phi_dir + i_phi, int8)
        end if
      end do

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

end module graph_module
