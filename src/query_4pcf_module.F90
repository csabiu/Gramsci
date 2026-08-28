module query_4pcf_module
  use kdtree2_precision_module
  use iso_fortran_env, only: int8, int64
  use config_module
  use graph_utils_module, only: find_dist
  implicit none

  real(kdkind), parameter :: PI = 3.141592653589793238462643383279502884197d0

  ! Tetrahedra whose direction-pixel triple product is (mathematically) zero
  ! have no defined chirality and contribute 0 to the parity-odd channel.
  ! Pixelized directions make exact-zero volumes common (repeated pixels,
  ! equatorial/coplanar pixel triples); deciding their sign from the floating-
  ! point residue would make the odd channel compiler-dependent noise.  The
  ! smallest GENUINELY nonzero pixel-triple volume shrinks with grid
  ! resolution -- measured 3.4e-3 at 4x16, 1.2e-4 at 8x32, 4.6e-6 at 16x64,
  ! 4.3e-8 at 23x91 and 5.2e-8 at the 128x256 int16 cap -- so the fixed
  ! tolerance stays at least ~40x below any legitimate configuration for every
  ! grid this build accepts, while sitting ~1e7 above fp noise (~1e-16).
  ! Under -exactparity it only ever catches genuinely coplanar quadruplets.
  real(kdkind), parameter :: VOL_DEGEN_TOL = 1.0d-9

  ! ---------------------------------------------------------------------------
  ! S4 symmetry tables for the 24 vertex permutations acting on 6 edges.
  !
  ! Edge numbering: 1=(1,2) 2=(1,3) 3=(1,4) 4=(2,3) 5=(2,4) 6=(3,4)
  !
  ! S4_EDGE_PERMS(:,k) gives the edge permutation for the k-th element of S4.
  ! Reading: the new 6-tuple is obtained by reading positions perm(:) of the
  ! original tuple, i.e. new(j) = old(perm(j)).
  !
  ! S4_PARITY(k) = +1 for even vertex permutations, -1 for odd.
  ! Even permutations preserve the sign of the signed volume (triple scalar
  ! product), odd permutations flip it.
  ! ---------------------------------------------------------------------------
  integer, parameter :: S4_EDGE_PERMS(6, 24) = reshape([ &
    1,2,3,4,5,6, &   !  1. ()           even
    1,4,5,2,3,6, &   !  2. (12)         odd
    4,2,6,1,5,3, &   !  3. (13)         odd
    5,6,3,4,1,2, &   !  4. (14)         odd
    2,1,3,4,6,5, &   !  5. (23)         odd
    3,2,1,6,5,4, &   !  6. (24)         odd
    1,3,2,5,4,6, &   !  7. (34)         odd
    1,5,4,3,2,6, &   !  8. (12)(34)     even
    6,2,4,3,5,1, &   !  9. (13)(24)     even
    6,5,3,4,2,1, &   ! 10. (14)(23)     even
    4,1,5,2,6,3, &   ! 11. (123)        even
    5,4,1,6,3,2, &   ! 12. (124)        even
    2,4,6,1,3,5, &   ! 13. (132)        even
    4,6,2,5,1,3, &   ! 14. (134)        even
    3,6,5,2,1,4, &   ! 15. (142)        even
    5,3,6,1,4,2, &   ! 16. (143)        even
    2,3,1,6,4,5, &   ! 17. (234)        even
    3,1,2,5,6,4, &   ! 18. (243)        even
    4,5,1,6,2,3, &   ! 19. (1234)       odd
    5,1,4,3,6,2, &   ! 20. (1243)       odd
    6,4,2,5,3,1, &   ! 21. (1324)       odd
    2,6,4,3,1,5, &   ! 22. (1342)       odd
    6,3,5,2,4,1, &   ! 23. (1423)       odd
    3,5,6,1,2,4  &   ! 24. (1432)       odd
    ], [6, 24])

  integer, parameter :: S4_PARITY(24) = [ &
     1, &   !  1. ()
    -1, &   !  2. (12)
    -1, &   !  3. (13)
    -1, &   !  4. (14)
    -1, &   !  5. (23)
    -1, &   !  6. (24)
    -1, &   !  7. (34)
     1, &   !  8. (12)(34)
     1, &   !  9. (13)(24)
     1, &   ! 10. (14)(23)
     1, &   ! 11. (123)
     1, &   ! 12. (124)
     1, &   ! 13. (132)
     1, &   ! 14. (134)
     1, &   ! 15. (142)
     1, &   ! 16. (143)
     1, &   ! 17. (234)
     1, &   ! 18. (243)
    -1, &   ! 19. (1234)
    -1, &   ! 20. (1243)
    -1, &   ! 21. (1324)
    -1, &   ! 22. (1342)
    -1, &   ! 23. (1423)
    -1  &   ! 24. (1432)
    ]

  ! ---------------------------------------------------------------------------
  ! Precomputed unit direction vectors for each direction pixel.
  ! Indexed by pixel index (1 to n_dir_pixels).
  ! Computed once at initialisation by init_direction_lookup().
  ! ---------------------------------------------------------------------------
  real(kdkind), allocatable :: dir_x(:), dir_y(:), dir_z(:)
  ! Position snapshot for -exactparity (filled before points is deallocated)
  real(kdkind), allocatable :: px(:), py(:), pz(:)

contains

  ! ---------------------------------------------------------------------------
  ! Initialise the direction pixel lookup table.
  ! Precomputes unit vectors at the centre of each (theta, phi) bin.
  ! ---------------------------------------------------------------------------
  subroutine init_direction_lookup()
    integer :: it, ip, idx
    real(kdkind) :: theta_c, phi_c

    allocate(dir_x(cfg%n_dir_pixels))
    allocate(dir_y(cfg%n_dir_pixels))
    allocate(dir_z(cfg%n_dir_pixels))

    do it = 1, cfg%n_theta_dir
      do ip = 1, cfg%n_phi_dir
        idx = (it - 1) * cfg%n_phi_dir + ip
        theta_c = (dble(it) - 0.5d0) * PI / dble(cfg%n_theta_dir)
        phi_c   = (dble(ip) - 0.5d0) * 2.0d0 * PI / dble(cfg%n_phi_dir) - PI
        dir_x(idx) = sin(theta_c) * cos(phi_c)
        dir_y(idx) = sin(theta_c) * sin(phi_c)
        dir_z(idx) = cos(theta_c)
      end do
    end do

    if (cfg%rank == 0) then
      print '(" 4PCF parity: direction pixels = ",i0," (n_theta=",i0, &
             &", n_phi=",i0,")")', &
             cfg%n_dir_pixels, cfg%n_theta_dir, cfg%n_phi_dir
    end if
  end subroutine init_direction_lookup

  ! ---------------------------------------------------------------------------
  ! Snapshot galaxy+random positions into 1D arrays for -exactparity.
  ! Must be called after graph construction, before points is deallocated.
  ! ---------------------------------------------------------------------------
  subroutine init_exact_parity_positions()
    integer :: n
    n = cfg%num_data + cfg%num_rand
    allocate(px(n), py(n), pz(n))
    px = points(1, 1:n)
    py = points(2, 1:n)
    pz = points(3, 1:n)
    if (cfg%rank == 0) &
      print *, '4PCF parity: exact signed volume from positions (-exactparity)'
  end subroutine init_exact_parity_positions

  ! ---------------------------------------------------------------------------
  ! Build the 6D lookup table bintable6(n,n,n,n,n,n) where n = cfg%nbins.
  ! ---------------------------------------------------------------------------
  subroutine create_4pcf_binlookup()
    integer :: b1, b2, b3, b4, b5, b6, n, k, n_configs_upper
    integer(int64) :: n64, upper64
    integer :: bins(6), perm_bins(6), canon(6)
    integer :: canon_parity, best_parity
    logical :: is_less

    n = cfg%nbins

    ! Burnside upper bound on the orbit count, in int64: n**6 overflows
    ! default integers already at nbins = 36 (36^6 > 2^31), which would turn
    ! the allocation size negative.  In practice memory runs out well before
    ! the int32 ceiling, but fail with an explanation rather than a crash.
    n64 = int(n, int64)
    upper64 = (n64**6 + 9_int64 * n64**4 + 14_int64 * n64**2) / 24_int64
    if (upper64 > int(huge(1), int64)) then
      print '("ERROR: -nbins ",i0," gives ~",i0," 4PCF configurations; reduce -nbins")', &
        n, upper64
      stop 1
    end if
    n_configs_upper = int(upper64)

    allocate(bintable6(n, n, n, n, n, n))
    bintable6 = 0
    cfg%n_configs_4pcf = 0

    ! Pre-allocate canon_bins_4pcf using the Burnside upper bound so we can
    ! record the canonical tuple for each new orbit as we encounter it.
    allocate(canon_bins_4pcf(6, n_configs_upper))
    allocate(orbit_mult_4pcf(n_configs_upper))
    orbit_mult_4pcf = 0

    do b1 = 1, n
    do b2 = 1, n
    do b3 = 1, n
    do b4 = 1, n
    do b5 = 1, n
    do b6 = 1, n
      if (bintable6(b1,b2,b3,b4,b5,b6) /= 0) cycle

      bins = [b1, b2, b3, b4, b5, b6]

      ! Find canonical (lexicographically smallest) permuted form
      canon = bins
      best_parity = S4_PARITY(1)  ! identity

      do k = 2, 24
        perm_bins = [ bins(S4_EDGE_PERMS(1,k)), bins(S4_EDGE_PERMS(2,k)), &
                      bins(S4_EDGE_PERMS(3,k)), bins(S4_EDGE_PERMS(4,k)), &
                      bins(S4_EDGE_PERMS(5,k)), bins(S4_EDGE_PERMS(6,k)) ]
        call lex_less(perm_bins, canon, is_less)
        if (is_less) then
          canon = perm_bins
          best_parity = S4_PARITY(k)
        end if
      end do

      ! Assign a new config index for this canonical form
      cfg%n_configs_4pcf = cfg%n_configs_4pcf + 1

      ! Record the canonical 6-tuple so write routines can iterate by config index
      canon_bins_4pcf(:, cfg%n_configs_4pcf) = canon

      ! Fill all 24 permutations of the original 6-tuple with this config index
      do k = 1, 24
        perm_bins = [ bins(S4_EDGE_PERMS(1,k)), bins(S4_EDGE_PERMS(2,k)), &
                      bins(S4_EDGE_PERMS(3,k)), bins(S4_EDGE_PERMS(4,k)), &
                      bins(S4_EDGE_PERMS(5,k)), bins(S4_EDGE_PERMS(6,k)) ]

        if (bintable6(perm_bins(1), perm_bins(2), perm_bins(3), &
                      perm_bins(4), perm_bins(5), perm_bins(6)) /= 0) cycle

        canon_parity = best_parity * S4_PARITY(k)

        bintable6(perm_bins(1), perm_bins(2), perm_bins(3), &
                  perm_bins(4), perm_bins(5), perm_bins(6)) = &
          canon_parity * cfg%n_configs_4pcf
        ! Count the distinct ordered 6-tuples in this orbit (each table
        ! entry is set exactly once); used by the analytic RRRR.
        orbit_mult_4pcf(cfg%n_configs_4pcf) = orbit_mult_4pcf(cfg%n_configs_4pcf) + 1
      end do

    end do
    end do
    end do
    end do
    end do
    end do

    if (cfg%rank == 0) then
      print *, '4PCF: number of unique configurations =', cfg%n_configs_4pcf
      print *, '4PCF: expected from Burnside formula  =', &
        (n**6 + 9*n**4 + 14*n**2) / 24
    end if

  end subroutine create_4pcf_binlookup

  ! ---------------------------------------------------------------------------
  ! Lexicographic comparison: is a(:) < b(:)?
  ! ---------------------------------------------------------------------------
  pure subroutine lex_less(a, b, result)
    integer, intent(in) :: a(6), b(6)
    logical, intent(out) :: result
    integer :: j
    result = .false.
    do j = 1, 6
      if (a(j) < b(j)) then
        result = .true.
        return
      else if (a(j) > b(j)) then
        return
      end if
    end do
  end subroutine lex_less

  ! ---------------------------------------------------------------------------
  ! Compute the isotropic 2PCF internally for disconnected 4PCF subtraction.
  !
  ! Accumulates pair counts into DD_2pcf(nbins) and RR_2pcf(nbins), then
  ! computes xi_2pcf(bin) = DD_2pcf(bin) / RR_2pcf(bin) (signed-weight
  ! numerator, so no "-1"; see note at the estimator below).
  ! Uses the same hub-based loop as query_graph_2pcf but with nmu=1 (isotropic)
  ! and dedicated arrays that do not interfere with N2/N3.
  ! ---------------------------------------------------------------------------
  subroutine compute_2pcf_for_4pcf(istart, iend)
    integer, intent(in) :: istart, iend
    integer :: i, k1, nn2, id1, bin_idx
    integer :: jr1, jr2
    integer(int8) :: ind1
    real(kdkind) :: wpair

    if (.not. allocated(DD_2pcf)) allocate(DD_2pcf(cfg%nbins))
    if (.not. allocated(RR_2pcf)) allocate(RR_2pcf(cfg%nbins))
    if (.not. allocated(xi_2pcf)) allocate(xi_2pcf(cfg%nbins))
    if (.not. allocated(xi2_2pcf)) allocate(xi2_2pcf(cfg%nbins))
    if (.not. allocated(xi4_2pcf)) allocate(xi4_2pcf(cfg%nbins))
    ! Jackknife touching sums of the pair counts; always allocated (dummy
    ! slice when off) because the OpenMP reduction clause names them.
    if (.not. allocated(DD_2pcf_jk)) allocate(DD_2pcf_jk(cfg%nbins, max(cfg%njk, 1)))
    if (.not. allocated(RR_2pcf_jk)) allocate(RR_2pcf_jk(cfg%nbins, max(cfg%njk, 1)))
    DD_2pcf = 0.0d0
    RR_2pcf = 0.0d0
    xi_2pcf = 0.0d0
    xi2_2pcf = 0.0d0
    xi4_2pcf = 0.0d0
    DD_2pcf_jk = 0.0d0
    RR_2pcf_jk = 0.0d0

    if (cfg%rank == 0) print *, 'Computing internal 2PCF for disconnected subtraction'

    !$OMP PARALLEL DO schedule(dynamic) &
    !$OMP& private(i, k1, nn2, id1, ind1, bin_idx, jr1, jr2, wpair) &
    !$OMP& shared(weights, output, buffer, cfg, region) &
    !$OMP& reduction(+:DD_2pcf, RR_2pcf) &
    !$OMP& reduction(+:DD_2pcf_jk, RR_2pcf_jk)
    do i = istart, iend
      if (buffer(i) == 1) cycle

      nn2 = output(i)%nn

      do k1 = 1, nn2
        ind1 = output(i)%dist(k1)
        id1 = output(i)%id(k1)
        bin_idx = int(ind1)
        wpair = weights(i) * weights(id1)

        ! All pairs (DD - 2DR + RR with signed weights)
        DD_2pcf(bin_idx) = DD_2pcf(bin_idx) + wpair

        ! Jackknife: each pair against each DISTINCT region it touches, so
        ! the per-realisation xi0 in the jackknifed disconnected term uses
        ! the same delete-one convention as the tuple counts.
        if (cfg%njk > 0) then
          jr1 = region(i)
          jr2 = region(id1)
          if (jr1 > 0) &
            DD_2pcf_jk(bin_idx, jr1) = DD_2pcf_jk(bin_idx, jr1) + wpair
          if (jr2 > 0 .and. jr2 /= jr1) &
            DD_2pcf_jk(bin_idx, jr2) = DD_2pcf_jk(bin_idx, jr2) + wpair
        end if

        ! RR pairs only (both points are randoms)
        if (i > cfg%num_data .and. id1 > cfg%num_data) then
          RR_2pcf(bin_idx) = RR_2pcf(bin_idx) + wpair
          if (cfg%njk > 0) then
            if (jr1 > 0) &
              RR_2pcf_jk(bin_idx, jr1) = RR_2pcf_jk(bin_idx, jr1) + wpair
            if (jr2 > 0 .and. jr2 /= jr1) &
              RR_2pcf_jk(bin_idx, jr2) = RR_2pcf_jk(bin_idx, jr2) + wpair
          end if
        end if
      end do
    end do
    !$OMP END PARALLEL DO

    ! Analytic mode: no random points exist, so RR comes from the periodic
    ! shell volumes and DD/RR - 1 is the natural estimator.  Otherwise
    ! xi = N2/RR with NO "-1": DD_2pcf is the signed all-pairs numerator
    ! (DD - 2DR + RR) under the data=+1/random=-1 weight convention, so it is
    ! already the correlation; it vanishes for an unclustered field.  This
    ! matches the standalone estimator write_2pcf_results.
    if (cfg%analytic) then
      call fill_rr_analytic()
    end if
    do bin_idx = 1, cfg%nbins
      if (RR_2pcf(bin_idx) /= 0.0d0) then
        if (cfg%analytic) then
          xi_2pcf(bin_idx) = DD_2pcf(bin_idx) / RR_2pcf(bin_idx) - 1.0d0
        else
          xi_2pcf(bin_idx) = DD_2pcf(bin_idx) / RR_2pcf(bin_idx)
        end if
      else
        xi_2pcf(bin_idx) = 0.0d0
      end if
    end do

    ! Anisotropic disconnected subtraction: xi_ell = (2*ell+1) * S_ell / RR,
    ! with S_ell the pair Legendre sums accumulated in create_graph.  Under
    ! the signed-weight convention S_ell is already the correlation-only
    ! part (the uniform term integrates to zero against L_ell); same in
    ! analytic mode, where no "-1" appears for ell > 0.
    if (cfg%disc_rsd) then
      do bin_idx = 1, cfg%nbins
        if (RR_2pcf(bin_idx) /= 0.0d0) then
          xi2_2pcf(bin_idx) = 5.0d0 * sum_pair_l2(bin_idx) / RR_2pcf(bin_idx)
          xi4_2pcf(bin_idx) = 9.0d0 * sum_pair_l4(bin_idx) / RR_2pcf(bin_idx)
        end if
      end do
    end if

    if (cfg%rank == 0) then
      print *, 'Internal 2PCF computed for disconnected subtraction:'
      do bin_idx = 1, cfg%nbins
        if (cfg%disc_rsd) then
          print '("  bin ",i3,": xi0 = ",e14.7,"  xi2 = ",e14.7,"  xi4 = ",e14.7)', &
            bin_idx, xi_2pcf(bin_idx), xi2_2pcf(bin_idx), xi4_2pcf(bin_idx)
        else
          print '("  bin ",i3,": xi = ",e14.7)', bin_idx, xi_2pcf(bin_idx)
        end if
      end do
    end if

  end subroutine compute_2pcf_for_4pcf

  ! Fill RR_2pcf with the analytic periodic-box shell counts.
  subroutine fill_rr_analytic()
    use analytic_randoms_module, only: analytic_rr
    call analytic_rr(RR_2pcf)
  end subroutine fill_rr_analytic

  ! ---------------------------------------------------------------------------
  ! All-configurations 4PCF query with parity decomposition.
  !
  ! Enumerates all tetrahedra (i, id1, id2, id3) using the hub-based approach.
  ! Uses graph-stored direction pixel indices for parity computation.
  !
  ! Result arrays: N4(config, 1) = total count (even channel)
  !                N4(config, 2) = parity-weighted count (odd channel)
  !                R4 = same but for RRRR (all-random) quadruplets
  ! ---------------------------------------------------------------------------
  ! Original binary-search implementation, kept for benchmarking.
  subroutine query_graph_4pcf_parity_bsearch(istart, iend)
    integer, intent(in) :: istart, iend
    integer :: i, k1, k2, k3, nn2, id1, id2, id3
    integer(int8) :: ind1, ind2, ind3, ind4, ind5, ind6
    integer :: raw_bin, config_idx, parity_flip, sign_V
    integer :: p1, p2, p3
    real(kdkind) :: w4, vol
    real(kdkind) :: u1x, u1y, u1z, u2x, u2y, u2z, u3x, u3y, u3z, rn
    logical :: exact_g

    if (cfg%njk > 0) then
      print *, 'ERROR: -njk is implemented for the merge-walk 4PCF only, '// &
               'not the binary-search variant.'
      stop
    end if
    exact_g = cfg%exact_parity
    if (cfg%rank == 0) print *, 'Performing 4PCF parity (binary-search)'
    if (cfg%rank == 0) print *, 'begin querying the graph'

    !$OMP PARALLEL DO schedule(dynamic) &
    !$OMP& private(i, k1, k2, k3, nn2, id1, id2, id3, &
    !$OMP&         ind1, ind2, ind3, ind4, ind5, ind6, &
    !$OMP&         raw_bin, config_idx, parity_flip, sign_V, &
    !$OMP&         p1, p2, p3, w4, vol, &
    !$OMP&         u1x, u1y, u1z, u2x, u2y, u2z, u3x, u3y, u3z, rn) &
    !$OMP& shared(weights, output, buffer, cfg, bintable6, dir_x, dir_y, dir_z, &
    !$OMP&        px, py, pz, exact_g) &
    !$OMP& reduction(+:N4, R4)
    do i = istart, iend
      if (buffer(i) == 1) cycle

      nn2 = output(i)%nn
      if (nn2 <= 2) cycle

      do k1 = 1, nn2
        ind1 = output(i)%dist(k1)
        id1 = output(i)%id(k1)

        do k2 = k1 + 1, nn2
          ind2 = output(i)%dist(k2)
          id2 = output(i)%id(k2)

          call find_dist(id1, id2, ind4)
          if (ind4 == 0) cycle

          do k3 = k2 + 1, nn2
            ind3 = output(i)%dist(k3)
            id3 = output(i)%id(k3)

            call find_dist(id1, id3, ind5)
            if (ind5 == 0) cycle

            call find_dist(id2, id3, ind6)
            if (ind6 == 0) cycle

            raw_bin = bintable6(ind1, ind2, ind3, ind4, ind5, ind6)
            config_idx = abs(raw_bin)
            parity_flip = sign(1, raw_bin)

            if (exact_g) then
              u1x = px(id1) - px(i); u1y = py(id1) - py(i); u1z = pz(id1) - pz(i)
              rn  = 1.0d0 / sqrt(u1x*u1x + u1y*u1y + u1z*u1z)
              u1x = u1x * rn; u1y = u1y * rn; u1z = u1z * rn
              u2x = px(id2) - px(i); u2y = py(id2) - py(i); u2z = pz(id2) - pz(i)
              rn  = 1.0d0 / sqrt(u2x*u2x + u2y*u2y + u2z*u2z)
              u2x = u2x * rn; u2y = u2y * rn; u2z = u2z * rn
              u3x = px(id3) - px(i); u3y = py(id3) - py(i); u3z = pz(id3) - pz(i)
              rn  = 1.0d0 / sqrt(u3x*u3x + u3y*u3y + u3z*u3z)
              u3x = u3x * rn; u3y = u3y * rn; u3z = u3z * rn
              vol = u1x * (u2y*u3z - u2z*u3y) &
                  + u1y * (u2z*u3x - u2x*u3z) &
                  + u1z * (u2x*u3y - u2y*u3x)
            else
              p1 = output(i)%phi(k1)
              p2 = output(i)%phi(k2)
              p3 = output(i)%phi(k3)
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

            w4 = weights(i) * weights(id1) * weights(id2) * weights(id3)

            N4(config_idx, 1) = N4(config_idx, 1) + w4
            N4(config_idx, 2) = N4(config_idx, 2) + (parity_flip * sign_V) * w4

            if (i > cfg%num_data .and. id1 > cfg%num_data .and. &
                id2 > cfg%num_data .and. id3 > cfg%num_data) then
              R4(config_idx, 1) = R4(config_idx, 1) + w4
              R4(config_idx, 2) = R4(config_idx, 2) + (parity_flip * sign_V) * w4
            end if

          end do
        end do
      end do
    end do
    !$OMP END PARALLEL DO

    call write_4pcf_results()
  end subroutine query_graph_4pcf_parity_bsearch
  ! (no jackknife write: the bsearch variant rejects -njk above)

  ! Merge-walk implementation (default).
  ! k2: 2-way walk — hub[k1+1..nn2] ∩ id1's list.
  ! k3: 3-way walk — hub[a+1..nn2] ∩ id1[b+1..] ∩ id2's list (entries > id2).
  ! Replaces three find_dist binary searches with pointer walks on sorted lists.
  subroutine query_graph_4pcf_parity(istart, iend)
    integer, intent(in) :: istart, iend
    integer :: i, k1, nn2, id1, id2
    integer :: nn_id1, nn_id2
    integer :: a, b, alpha, beta, gamma
    integer :: ha, hb, hc
    integer :: jrs(4), jr, q, p
    integer(int8) :: ind1, ind2, ind3, ind4, ind5, ind6
    integer :: raw_bin, config_idx, parity_flip, sign_V
    integer :: p1, p2, p3
    real(kdkind) :: w4, vol, sv4
    real(kdkind) :: u1x, u1y, u1z, u2x, u2y, u2z, u3x, u3y, u3z, rn
    logical :: exact_g, dup, isr

    if (.not. cfg%half_graph) then
      print *, 'ERROR: merge-walk 4PCF parity requires half_graph=.true.'
      stop 1
    end if
    exact_g = cfg%exact_parity
    if (cfg%rank == 0) print *, 'Performing 4PCF parity (all configurations, merge-walk)'
    if (cfg%rank == 0) print *, 'begin querying the graph'

    !$OMP PARALLEL DO schedule(dynamic) &
    !$OMP& private(i, k1, nn2, id1, id2, nn_id1, nn_id2, &
    !$OMP&         a, b, alpha, beta, gamma, ha, hb, hc, &
    !$OMP&         ind1, ind2, ind3, ind4, ind5, ind6, &
    !$OMP&         raw_bin, config_idx, parity_flip, sign_V, &
    !$OMP&         p1, p2, p3, w4, vol, sv4, &
    !$OMP&         jrs, jr, q, p, dup, isr, &
    !$OMP&         u1x, u1y, u1z, u2x, u2y, u2z, u3x, u3y, u3z, rn) &
    !$OMP& shared(weights, output, buffer, cfg, bintable6, dir_x, dir_y, dir_z, &
    !$OMP&        px, py, pz, exact_g, region, N4jk, R4jk) &
    !$OMP& reduction(+:N4, R4)
    do i = istart, iend
      if (buffer(i) == 1) cycle
      nn2 = output(i)%nn
      if (nn2 <= 2) cycle

      do k1 = 1, nn2
        ind1   = output(i)%dist(k1)
        id1    = output(i)%id(k1)
        nn_id1 = output(id1)%nn
        if (exact_g) then
          u1x = px(id1) - px(i)
          u1y = py(id1) - py(i)
          u1z = pz(id1) - pz(i)
          rn  = 1.0d0 / sqrt(u1x*u1x + u1y*u1y + u1z*u1z)
          u1x = u1x * rn; u1y = u1y * rn; u1z = u1z * rn
        end if

        ! k2: 2-way merge walk
        a = k1 + 1
        b = 1
        do while (a <= nn2 .and. b <= nn_id1)
          if (output(i)%id(a) == output(id1)%id(b)) then
            ind2   = output(i)%dist(a)
            ind4   = output(id1)%dist(b)
            id2    = output(i)%id(a)
            nn_id2 = output(id2)%nn
            if (exact_g) then
              u2x = px(id2) - px(i)
              u2y = py(id2) - py(i)
              u2z = pz(id2) - pz(i)
              rn  = 1.0d0 / sqrt(u2x*u2x + u2y*u2y + u2z*u2z)
              u2x = u2x * rn; u2y = u2y * rn; u2z = u2z * rn
            end if

            ! k3: 3-way merge walk
            alpha = a + 1
            beta  = b + 1
            gamma = 1
            ! Fortran .and. does not short-circuit: guard the id() access
            do while (gamma <= nn_id2)
              if (output(id2)%id(gamma) > id2) exit
              gamma = gamma + 1
            end do

            do while (alpha <= nn2 .and. beta <= nn_id1 .and. gamma <= nn_id2)
              ha = output(i)%id(alpha)
              hb = output(id1)%id(beta)
              hc = output(id2)%id(gamma)
              if (ha == hb .and. ha == hc) then
                ind3 = output(i)%dist(alpha)
                ind5 = output(id1)%dist(beta)
                ind6 = output(id2)%dist(gamma)
                raw_bin     = bintable6(ind1, ind2, ind3, ind4, ind5, ind6)
                config_idx  = abs(raw_bin)
                parity_flip = sign(1, raw_bin)

                if (exact_g) then
                  u3x = px(ha) - px(i)
                  u3y = py(ha) - py(i)
                  u3z = pz(ha) - pz(i)
                  rn  = 1.0d0 / sqrt(u3x*u3x + u3y*u3y + u3z*u3z)
                  u3x = u3x * rn; u3y = u3y * rn; u3z = u3z * rn
                  vol = u1x * (u2y*u3z - u2z*u3y) &
                      + u1y * (u2z*u3x - u2x*u3z) &
                      + u1z * (u2x*u3y - u2y*u3x)
                else
                  p1 = output(i)%phi(k1)
                  p2 = output(i)%phi(a)
                  p3 = output(i)%phi(alpha)
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

                w4 = weights(i) * weights(id1) * weights(id2) * weights(ha)

                N4(config_idx, 1) = N4(config_idx, 1) + w4
                N4(config_idx, 2) = N4(config_idx, 2) + (parity_flip * sign_V) * w4

                isr = i > cfg%num_data .and. id1 > cfg%num_data .and. &
                      id2 > cfg%num_data .and. ha > cfg%num_data
                if (isr) then
                  R4(config_idx, 1) = R4(config_idx, 1) + w4
                  R4(config_idx, 2) = R4(config_idx, 2) + (parity_flip * sign_V) * w4
                end if

                ! Jackknife: accumulate this quadruplet against each DISTINCT
                ! region it touches (dedup as in the 3PCF).  ATOMIC, not a
                ! reduction: see allocate_4pcf_jk.
                if (cfg%njk > 0) then
                  jrs(1) = region(i)
                  jrs(2) = region(id1)
                  jrs(3) = region(id2)
                  jrs(4) = region(ha)
                  sv4 = (parity_flip * sign_V) * w4
                  do q = 1, 4
                    jr = jrs(q)
                    if (jr <= 0) cycle
                    dup = .false.
                    do p = 1, q - 1
                      if (jrs(p) == jr) dup = .true.
                    end do
                    if (dup) cycle
                    !$OMP ATOMIC UPDATE
                    N4jk(config_idx, 1, jr) = N4jk(config_idx, 1, jr) + w4
                    !$OMP ATOMIC UPDATE
                    N4jk(config_idx, 2, jr) = N4jk(config_idx, 2, jr) + sv4
                    if (isr) then
                      !$OMP ATOMIC UPDATE
                      R4jk(config_idx, 1, jr) = R4jk(config_idx, 1, jr) + w4
                      !$OMP ATOMIC UPDATE
                      R4jk(config_idx, 2, jr) = R4jk(config_idx, 2, jr) + sv4
                    end if
                  end do
                end if

                alpha = alpha + 1
                beta  = beta  + 1
                gamma = gamma + 1
              else
                if (ha <= hb .and. ha <= hc) alpha = alpha + 1
                if (hb <= ha .and. hb <= hc) beta  = beta  + 1
                if (hc <= ha .and. hc <= hb) gamma = gamma + 1
              end if
            end do

            a = a + 1
            b = b + 1
          else if (output(i)%id(a) < output(id1)%id(b)) then
            a = a + 1
          else
            b = b + 1
          end if
        end do
      end do
    end do
    !$OMP END PARALLEL DO

    call write_4pcf_results()
    call write_4pcf_jackknife(.true.)
  end subroutine query_graph_4pcf_parity

  ! ---------------------------------------------------------------------------
  ! All-configurations 4PCF query WITHOUT parity decomposition.
  !
  ! Same enumeration as the parity version but does not use direction pixels.
  ! Accumulates into a single count channel: N4(config, 1) and R4(config, 1).
  ! ---------------------------------------------------------------------------
  ! Original binary-search implementation, kept for benchmarking.
  subroutine query_graph_4pcf_bsearch(istart, iend)
    integer, intent(in) :: istart, iend
    integer :: i, k1, k2, k3, nn2, id1, id2, id3
    integer(int8) :: ind1, ind2, ind3, ind4, ind5, ind6
    integer :: config_idx
    real(kdkind) :: w4

    if (cfg%njk > 0) then
      print *, 'ERROR: -njk is implemented for the merge-walk 4PCF only, '// &
               'not the binary-search variant.'
      stop
    end if
    if (cfg%rank == 0) print *, 'Performing 4PCF (binary-search)'
    if (cfg%rank == 0) print *, 'begin querying the graph'

    !$OMP PARALLEL DO schedule(dynamic) &
    !$OMP& private(i, k1, k2, k3, nn2, id1, id2, id3, &
    !$OMP&         ind1, ind2, ind3, ind4, ind5, ind6, &
    !$OMP&         config_idx, w4) &
    !$OMP& shared(weights, output, buffer, cfg, bintable6) &
    !$OMP& reduction(+:N4, R4)
    do i = istart, iend
      if (buffer(i) == 1) cycle

      nn2 = output(i)%nn
      if (nn2 <= 2) cycle

      do k1 = 1, nn2
        ind1 = output(i)%dist(k1)
        id1 = output(i)%id(k1)

        do k2 = k1 + 1, nn2
          ind2 = output(i)%dist(k2)
          id2 = output(i)%id(k2)

          call find_dist(id1, id2, ind4)
          if (ind4 == 0) cycle

          do k3 = k2 + 1, nn2
            ind3 = output(i)%dist(k3)
            id3 = output(i)%id(k3)

            call find_dist(id1, id3, ind5)
            if (ind5 == 0) cycle

            call find_dist(id2, id3, ind6)
            if (ind6 == 0) cycle

            config_idx = abs(bintable6(ind1, ind2, ind3, ind4, ind5, ind6))
            if (config_idx == 0) cycle

            w4 = weights(i) * weights(id1) * weights(id2) * weights(id3)

            N4(config_idx, 1) = N4(config_idx, 1) + w4

            if (i > cfg%num_data .and. id1 > cfg%num_data .and. &
                id2 > cfg%num_data .and. id3 > cfg%num_data) then
              R4(config_idx, 1) = R4(config_idx, 1) + w4
            end if

          end do
        end do
      end do
    end do
    !$OMP END PARALLEL DO

    call write_4pcf_results_noparity()
  end subroutine query_graph_4pcf_bsearch

  ! Merge-walk implementation (default).
  ! k2: 2-way walk — hub[k1+1..nn2] ∩ id1's list.
  ! k3: 3-way walk — hub[a+1..nn2] ∩ id1[b+1..] ∩ id2's list (entries > id2).
  ! Replaces three find_dist binary searches with pointer walks on sorted lists.
  subroutine query_graph_4pcf(istart, iend)
    integer, intent(in) :: istart, iend
    integer :: i, k1, nn2, id1, id2
    integer :: nn_id1, nn_id2
    integer :: a, b, alpha, beta, gamma
    integer :: ha, hb, hc, config_idx
    integer :: jrs(4), jr, q, p
    integer(int8) :: ind1, ind2, ind3, ind4, ind5, ind6
    real(kdkind) :: w4
    logical :: dup, isr

    if (.not. cfg%half_graph) then
      print *, 'ERROR: merge-walk 4PCF requires half_graph=.true.'
      stop 1
    end if
    if (cfg%rank == 0) print *, 'Performing 4PCF (all configurations, merge-walk)'
    if (cfg%rank == 0) print *, 'begin querying the graph'

    !$OMP PARALLEL DO schedule(dynamic) &
    !$OMP& private(i, k1, nn2, id1, id2, nn_id1, nn_id2, &
    !$OMP&         a, b, alpha, beta, gamma, ha, hb, hc, config_idx, &
    !$OMP&         jrs, jr, q, p, dup, isr, &
    !$OMP&         ind1, ind2, ind3, ind4, ind5, ind6, w4) &
    !$OMP& shared(weights, output, buffer, cfg, bintable6, region, N4jk, R4jk) &
    !$OMP& reduction(+:N4, R4)
    do i = istart, iend
      if (buffer(i) == 1) cycle
      nn2 = output(i)%nn
      if (nn2 <= 2) cycle

      do k1 = 1, nn2
        ind1   = output(i)%dist(k1)
        id1    = output(i)%id(k1)
        nn_id1 = output(id1)%nn

        ! k2: 2-way merge walk
        a = k1 + 1
        b = 1
        do while (a <= nn2 .and. b <= nn_id1)
          if (output(i)%id(a) == output(id1)%id(b)) then
            ind2   = output(i)%dist(a)
            ind4   = output(id1)%dist(b)
            id2    = output(i)%id(a)
            nn_id2 = output(id2)%nn

            ! k3: 3-way merge walk
            alpha = a + 1
            beta  = b + 1
            gamma = 1
            ! Fortran .and. does not short-circuit: guard the id() access
            do while (gamma <= nn_id2)
              if (output(id2)%id(gamma) > id2) exit
              gamma = gamma + 1
            end do

            do while (alpha <= nn2 .and. beta <= nn_id1 .and. gamma <= nn_id2)
              ha = output(i)%id(alpha)
              hb = output(id1)%id(beta)
              hc = output(id2)%id(gamma)
              if (ha == hb .and. ha == hc) then
                ind3 = output(i)%dist(alpha)
                ind5 = output(id1)%dist(beta)
                ind6 = output(id2)%dist(gamma)
                config_idx = abs(bintable6(ind1, ind2, ind3, ind4, ind5, ind6))
                if (config_idx /= 0) then
                  w4 = weights(i) * weights(id1) * weights(id2) * weights(ha)
                  N4(config_idx, 1) = N4(config_idx, 1) + w4
                  isr = i > cfg%num_data .and. id1 > cfg%num_data .and. &
                        id2 > cfg%num_data .and. ha > cfg%num_data
                  if (isr) then
                    R4(config_idx, 1) = R4(config_idx, 1) + w4
                  end if
                  ! Jackknife: each quadruplet against each DISTINCT region
                  ! it touches; ATOMIC, not a reduction (allocate_4pcf_jk).
                  if (cfg%njk > 0) then
                    jrs(1) = region(i)
                    jrs(2) = region(id1)
                    jrs(3) = region(id2)
                    jrs(4) = region(ha)
                    do q = 1, 4
                      jr = jrs(q)
                      if (jr <= 0) cycle
                      dup = .false.
                      do p = 1, q - 1
                        if (jrs(p) == jr) dup = .true.
                      end do
                      if (dup) cycle
                      !$OMP ATOMIC UPDATE
                      N4jk(config_idx, 1, jr) = N4jk(config_idx, 1, jr) + w4
                      if (isr) then
                        !$OMP ATOMIC UPDATE
                        R4jk(config_idx, 1, jr) = R4jk(config_idx, 1, jr) + w4
                      end if
                    end do
                  end if
                end if
                alpha = alpha + 1
                beta  = beta  + 1
                gamma = gamma + 1
              else
                if (ha <= hb .and. ha <= hc) alpha = alpha + 1
                if (hb <= ha .and. hb <= hc) beta  = beta  + 1
                if (hc <= ha .and. hc <= hb) gamma = gamma + 1
              end if
            end do

            a = a + 1
            b = b + 1
          else if (output(i)%id(a) < output(id1)%id(b)) then
            a = a + 1
          else
            b = b + 1
          end if
        end do
      end do
    end do
    !$OMP END PARALLEL DO

    call write_4pcf_results_noparity()
    call write_4pcf_jackknife(.false.)
  end subroutine query_graph_4pcf

  ! ---------------------------------------------------------------------------
  ! Disconnected (Gaussian) 4PCF for one canonical bin 6-tuple: the sum over
  ! the 3 complementary edge pairings of the orientation-averaged product
  ! <xi(r_a, mu_a) xi(r_b, mu_b)>.  Because the two opposite edges are rigidly
  ! attached to the same tetrahedron, their line-of-sight angles co-vary;
  ! by the Legendre addition theorem the average is
  !   xi0*xi0 + xi2*xi2 * L2(ct)/5 + xi4*xi4 * L4(ct)/9
  ! with ct the cosine of the inter-edge angle, fixed by the six side
  ! lengths (evaluated at bin centres):  e12.e34 = (r14^2+r23^2-r13^2-r24^2)/2.
  ! xi2/xi4 are zero unless cfg%disc_rsd accumulated pair multipoles, in
  ! which case this reduces to the isotropic product subtraction.
  ! ---------------------------------------------------------------------------
  ! xi0 is passed explicitly so the jackknife write can evaluate the term
  ! with each realisation's delete-one xi0 (the full-sample callers pass
  ! xi_2pcf).  xi2/xi4 stay at their module (full-sample) values.
  function zeta_disc_config(b, xi0) result(disc)
    integer, intent(in) :: b(6)
    real(kdkind), intent(in) :: xi0(:)
    real(kdkind) :: disc
    real(kdkind) :: rc(6), rc2(6), ct, l2ct, l4ct
    integer :: k, ea, eb
    ! Complementary edge pairings: (1,6)=(12,34), (2,5)=(13,24), (3,4)=(14,23)
    integer, parameter :: PAIR_A(3) = [1, 2, 3]
    integer, parameter :: PAIR_B(3) = [6, 5, 4]

    do k = 1, 6
      rc(k)  = 0.5d0 * (radial_bins(b(k)) + radial_bins(b(k) + 1))
      rc2(k) = rc(k) * rc(k)
    end do

    disc = 0.0d0
    do k = 1, 3
      select case (k)
      case (1)   ! edges (1,2) & (3,4)
        ct = (rc2(3) + rc2(4) - rc2(2) - rc2(5)) / (2.0d0 * rc(1) * rc(6))
      case (2)   ! edges (1,3) & (2,4)
        ct = (rc2(3) + rc2(4) - rc2(1) - rc2(6)) / (2.0d0 * rc(2) * rc(5))
      case (3)   ! edges (1,4) & (2,3)
        ct = (rc2(2) + rc2(5) - rc2(1) - rc2(6)) / (2.0d0 * rc(3) * rc(4))
      end select
      ! Bin centres of a degenerate config can slightly violate |ct| <= 1
      ct = max(-1.0d0, min(1.0d0, ct))
      l2ct = 1.5d0 * ct * ct - 0.5d0
      l4ct = 4.375d0 * ct**4 - 3.75d0 * ct * ct + 0.375d0

      ea = b(PAIR_A(k))
      eb = b(PAIR_B(k))
      disc = disc + xi0(ea) * xi0(eb) &
           + xi2_2pcf(ea) * xi2_2pcf(eb) * l2ct / 5.0d0 &
           + xi4_2pcf(ea) * xi4_2pcf(eb) * l4ct / 9.0d0
    end do
  end function zeta_disc_config

  ! ---------------------------------------------------------------------------
  ! Write 4PCF parity results to output file.
  ! ---------------------------------------------------------------------------
  subroutine write_4pcf_results()
    use analytic_randoms_module, only: analytic_rrrr
    integer :: config_idx, unit_num
    integer :: b1, b2, b3, b4, b5, b6
    real(kdkind) :: zeta_even, zeta_odd, zeta_disc, zeta_conn_even, zeta_conn_odd

    if (cfg%rank /= cfg%master) return

    ! Analytic mode: RRRR from the periodic-box tetrahedron kernel.  The
    ! parity-odd random count vanishes exactly (mirror quadruplets cancel).
    if (cfg%analytic) then
      call analytic_rrrr(R4(:, 1))
      R4(:, 2) = 0.0d0
    end if

    unit_num = 40
    open(unit_num, file=trim(mode_output_file('4pcfp')), status='unknown')
    call write_provenance(unit_num)
    write(unit_num, '(a)') '# r12min r12max r13min r13max r14min r14max ' // &
                           'r23min r23max r24min r24max r34min r34max ' // &
                           'NNNN RRRR zeta_even NNNN_odd RRRR_odd zeta_odd ' // &
                           'zeta_disc zeta_conn_even zeta_conn_odd'

    ! Iterate over all unique configurations by index; use the stored canonical
    ! bin tuple for edge labels.  This covers every orbit regardless of whether
    ! its canonical representative is a non-decreasing 6-tuple.
    do config_idx = 1, cfg%n_configs_4pcf
      b1 = canon_bins_4pcf(1, config_idx)
      b2 = canon_bins_4pcf(2, config_idx)
      b3 = canon_bins_4pcf(3, config_idx)
      b4 = canon_bins_4pcf(4, config_idx)
      b5 = canon_bins_4pcf(5, config_idx)
      b6 = canon_bins_4pcf(6, config_idx)

      if (R4(config_idx, 1) /= 0.0d0) then
        zeta_even = N4(config_idx, 1) / R4(config_idx, 1)
      else
        zeta_even = 0.0d0
      end if

      if (R4(config_idx, 1) /= 0.0d0) then
        zeta_odd = N4(config_idx, 2) / R4(config_idx, 1)
      else
        zeta_odd = 0.0d0
      end if

      ! Analytic mode: NNNN is pure DDDD, so NNNN/RRRR = 1 + sum_6 xi(edge)
      ! + sum_4 zeta3(face) + [xi*xi pairings + connected 4PCF].  Subtract
      ! the lower-order terms so zeta_even keeps the same meaning as in the
      ! signed-weight estimator.  The parity-weighted numerator needs no
      ! subtraction: every lower-order term is parity-even and cancels.
      if (cfg%analytic .and. R4(config_idx, 1) /= 0.0d0) then
        zeta_even = zeta_even - 1.0d0 &
          - xi_2pcf(b1) - xi_2pcf(b2) - xi_2pcf(b3) &
          - xi_2pcf(b4) - xi_2pcf(b5) - xi_2pcf(b6) &
          - zeta3_internal(bintable(b1, b2, b4, 1)) &
          - zeta3_internal(bintable(b1, b3, b5, 1)) &
          - zeta3_internal(bintable(b2, b3, b6, 1)) &
          - zeta3_internal(bintable(b4, b5, b6, 1))
      end if

      ! Disconnected 4PCF: 3 complementary edge pairings (parity-even),
      ! including the redshift-space multipole covariance when available.
      zeta_disc = zeta_disc_config([b1, b2, b3, b4, b5, b6], xi_2pcf)

      ! Connected even channel: subtract disconnected (which is parity-even).
      ! zeta_even = N4/R4 is already the total 4PCF (no "-1"; see no-parity
      ! write routine).
      zeta_conn_even = zeta_even - zeta_disc

      ! Connected odd channel: disconnected contribution is zero for odd parity
      zeta_conn_odd = zeta_odd

      write(unit_num, '(12(e14.7,1x), 9(e14.7,1x))') &
        radial_bins(b1), radial_bins(b1+1), &
        radial_bins(b2), radial_bins(b2+1), &
        radial_bins(b3), radial_bins(b3+1), &
        radial_bins(b4), radial_bins(b4+1), &
        radial_bins(b5), radial_bins(b5+1), &
        radial_bins(b6), radial_bins(b6+1), &
        N4(config_idx, 1), R4(config_idx, 1), zeta_even, &
        N4(config_idx, 2), R4(config_idx, 2), zeta_odd, &
        zeta_disc, zeta_conn_even, zeta_conn_odd
    end do
    close(unit_num)

    if (cfg%rank == 0) print *, '4PCF parity results written to ', trim(mode_output_file('4pcfp'))
  end subroutine write_4pcf_results

  ! ---------------------------------------------------------------------------
  ! Write 4PCF results without parity channels.
  ! ---------------------------------------------------------------------------
  subroutine write_4pcf_results_noparity()
    use analytic_randoms_module, only: analytic_rrrr
    integer :: config_idx, unit_num
    integer :: b1, b2, b3, b4, b5, b6
    real(kdkind) :: zeta, zeta_disc, zeta_conn

    if (cfg%rank /= cfg%master) return

    ! Analytic mode: RRRR from the periodic-box tetrahedron kernel.
    if (cfg%analytic) then
      call analytic_rrrr(R4(:, 1))
    end if

    unit_num = 40
    open(unit_num, file=trim(mode_output_file('4pcf')), status='unknown')
    call write_provenance(unit_num)
    write(unit_num, '(a)') '# r12min r12max r13min r13max r14min r14max ' // &
                           'r23min r23max r24min r24max r34min r34max ' // &
                           'NNNN RRRR zeta zeta_disc zeta_conn'

    ! Iterate over all unique configurations by index; use the stored canonical
    ! bin tuple for edge labels.  This covers every orbit regardless of whether
    ! its canonical representative is a non-decreasing 6-tuple.
    do config_idx = 1, cfg%n_configs_4pcf
      b1 = canon_bins_4pcf(1, config_idx)
      b2 = canon_bins_4pcf(2, config_idx)
      b3 = canon_bins_4pcf(3, config_idx)
      b4 = canon_bins_4pcf(4, config_idx)
      b5 = canon_bins_4pcf(5, config_idx)
      b6 = canon_bins_4pcf(6, config_idx)

      if (R4(config_idx, 1) /= 0.0d0) then
        zeta = N4(config_idx, 1) / R4(config_idx, 1)
      else
        zeta = 0.0d0
      end if

      ! Analytic mode: NNNN is pure DDDD, so NNNN/RRRR = 1 + sum_6 xi(edge)
      ! + sum_4 zeta3(face) + [xi*xi pairings + connected 4PCF].  Subtract
      ! the lower-order terms so zeta keeps the same meaning as in the
      ! signed-weight estimator.
      if (cfg%analytic .and. R4(config_idx, 1) /= 0.0d0) then
        zeta = zeta - 1.0d0 &
          - xi_2pcf(b1) - xi_2pcf(b2) - xi_2pcf(b3) &
          - xi_2pcf(b4) - xi_2pcf(b5) - xi_2pcf(b6) &
          - zeta3_internal(bintable(b1, b2, b4, 1)) &
          - zeta3_internal(bintable(b1, b3, b5, 1)) &
          - zeta3_internal(bintable(b2, b3, b6, 1)) &
          - zeta3_internal(bintable(b4, b5, b6, 1))
      end if

      ! Disconnected 4PCF: 3 complementary edge pairings
      ! {(1,2),(3,4)} = {b1,b6}, {(1,3),(2,4)} = {b2,b5}, {(1,4),(2,3)} = {b3,b4}
      ! Invariant under S4 vertex permutations; includes the redshift-space
      ! multipole covariance terms when available.
      zeta_disc = zeta_disc_config([b1, b2, b3, b4, b5, b6], xi_2pcf)

      ! Connected 4PCF: total minus disconnected.  zeta = N4/R4 with signed
      ! weights is already the total 4PCF (it vanishes for an unclustered
      ! field), so no "-1" is applied.
      zeta_conn = zeta - zeta_disc

      write(unit_num, '(12(e14.7,1x), 5(e14.7,1x))') &
        radial_bins(b1), radial_bins(b1+1), &
        radial_bins(b2), radial_bins(b2+1), &
        radial_bins(b3), radial_bins(b3+1), &
        radial_bins(b4), radial_bins(b4+1), &
        radial_bins(b5), radial_bins(b5+1), &
        radial_bins(b6), radial_bins(b6+1), &
        N4(config_idx, 1), R4(config_idx, 1), zeta, zeta_disc, zeta_conn
    end do
    close(unit_num)

    if (cfg%rank == 0) print *, '4PCF results written to ', trim(mode_output_file('4pcf'))
  end subroutine write_4pcf_results_noparity

  ! ---------------------------------------------------------------------------
  ! Delete-one jackknife realisations and error summary for the 4PCF
  ! (parity = .true. adds the odd channel).  Realisation m omits every
  ! quadruplet touching region m:  N_m = N_total - N_touching(m).  Weight
  ! renormalisation is deliberately NOT applied (the deleted-region factor
  ! cancels in the ratio; see the 3PCF jackknife note).
  !
  ! The disconnected term is rebuilt per realisation from the jackknifed
  ! internal 2PCF xi0 (DD_2pcf_jk / RR_2pcf_jk); the RSD multipoles xi2/xi4
  ! (accumulated at graph build) are kept at their full-sample values --
  ! they are a subdominant part of zeta_disc.  The parity-odd channel needs
  ! no disconnected subtraction (all lower-order terms are parity-even), so
  ! zeta_conn_odd_m = zeta_odd_m identically.
  !
  ! Only the CPU merge-walk kernels fill N4jk/R4jk; the GPU drivers route
  ! 4PCF queries to the CPU when -njk is active.
  ! ---------------------------------------------------------------------------
  subroutine write_4pcf_jackknife(parity)
    logical, intent(in) :: parity
    integer :: config_idx, m, bb, unit_num, unit_err, b(6)
    real(kdkind) :: ddm, rrm, n4m, r4m, n4om, disc
    real(kdkind) :: zeta_m(cfg%njk), conn_m(cfg%njk), odd_m(cfg%njk)
    real(kdkind) :: z_mean, z_sigma, c_mean, c_sigma, o_mean, o_sigma
    real(kdkind), allocatable :: xi0_all(:,:)
    character(len=5) :: mode

    if (cfg%rank /= cfg%master) return
    if (cfg%njk <= 0) return

    mode = '4pcf'
    if (parity) mode = '4pcfp'

    ! Per-realisation delete-one xi0 for the disconnected term.
    allocate(xi0_all(cfg%nbins, cfg%njk))
    do m = 1, cfg%njk
      do bb = 1, cfg%nbins
        ddm = DD_2pcf(bb) - DD_2pcf_jk(bb, m)
        rrm = RR_2pcf(bb) - RR_2pcf_jk(bb, m)
        if (rrm /= 0.0d0) then
          xi0_all(bb, m) = ddm / rrm
        else
          xi0_all(bb, m) = 0.0d0
        end if
      end do
    end do

    unit_num = 41
    unit_err = 43
    open(unit_num, file=trim(mode_output_file(trim(mode)))//'.jk', status='unknown')
    call write_provenance(unit_num)
    write(unit_num, '(a,i0)') '# delete-one jackknife realisations, njk = ', cfg%njk
    write(unit_num, '(a)') '# zeta_disc uses each realisation''s delete-one xi0;'
    write(unit_num, '(a)') '#   the RSD multipoles xi2/xi4 are full-sample.'
    if (parity) then
      write(unit_num, '(a)') '# r12min r12max r13min r13max r14min r14max ' // &
        'r23min r23max r24min r24max r34min r34max ' // &
        'ireal NNNN RRRR zeta_even NNNN_odd zeta_odd ' // &
        'zeta_disc zeta_conn_even zeta_conn_odd'
    else
      write(unit_num, '(a)') '# r12min r12max r13min r13max r14min r14max ' // &
        'r23min r23max r24min r24max r34min r34max ' // &
        'ireal NNNN RRRR zeta zeta_disc zeta_conn'
    end if
    open(unit_err, file=trim(mode_output_file(trim(mode)))//'.jkerr', status='unknown')
    call write_provenance(unit_err)
    write(unit_err, '(a,i0)') '# delete-one jackknife error, njk = ', cfg%njk
    if (parity) then
      write(unit_err, '(a)') '# r12min r12max r13min r13max r14min r14max ' // &
        'r23min r23max r24min r24max r34min r34max ' // &
        'zeta_even_mean_jk zeta_even_sigma_jk zeta_odd_mean_jk zeta_odd_sigma_jk ' // &
        'zeta_conn_even_mean_jk zeta_conn_even_sigma_jk'
    else
      write(unit_err, '(a)') '# r12min r12max r13min r13max r14min r14max ' // &
        'r23min r23max r24min r24max r34min r34max ' // &
        'zeta_mean_jk zeta_sigma_jk zeta_conn_mean_jk zeta_conn_sigma_jk'
    end if

    do config_idx = 1, cfg%n_configs_4pcf
      b = canon_bins_4pcf(:, config_idx)
      do m = 1, cfg%njk
        n4m = N4(config_idx, 1) - N4jk(config_idx, 1, m)
        r4m = R4(config_idx, 1) - R4jk(config_idx, 1, m)
        if (r4m /= 0.0d0) then
          zeta_m(m) = n4m / r4m
        else
          zeta_m(m) = 0.0d0
        end if
        disc = zeta_disc_config(b, xi0_all(:, m))
        conn_m(m) = zeta_m(m) - disc
        if (parity) then
          n4om = N4(config_idx, 2) - N4jk(config_idx, 2, m)
          if (r4m /= 0.0d0) then
            odd_m(m) = n4om / r4m
          else
            odd_m(m) = 0.0d0
          end if
          write(unit_num, '(12(e14.7,1x),i6,1x,8(e14.7,1x))') &
            radial_bins(b(1)), radial_bins(b(1)+1), &
            radial_bins(b(2)), radial_bins(b(2)+1), &
            radial_bins(b(3)), radial_bins(b(3)+1), &
            radial_bins(b(4)), radial_bins(b(4)+1), &
            radial_bins(b(5)), radial_bins(b(5)+1), &
            radial_bins(b(6)), radial_bins(b(6)+1), &
            m, n4m, r4m, zeta_m(m), n4om, odd_m(m), disc, conn_m(m), odd_m(m)
        else
          write(unit_num, '(12(e14.7,1x),i6,1x,5(e14.7,1x))') &
            radial_bins(b(1)), radial_bins(b(1)+1), &
            radial_bins(b(2)), radial_bins(b(2)+1), &
            radial_bins(b(3)), radial_bins(b(3)+1), &
            radial_bins(b(4)), radial_bins(b(4)+1), &
            radial_bins(b(5)), radial_bins(b(5)+1), &
            radial_bins(b(6)), radial_bins(b(6)+1), &
            m, n4m, r4m, zeta_m(m), disc, conn_m(m)
        end if
      end do
      call jk_mean_sigma(zeta_m, z_mean, z_sigma)
      call jk_mean_sigma(conn_m, c_mean, c_sigma)
      if (parity) then
        call jk_mean_sigma(odd_m, o_mean, o_sigma)
        write(unit_err, '(12(e14.7,1x),6(e14.7,1x))') &
          radial_bins(b(1)), radial_bins(b(1)+1), &
          radial_bins(b(2)), radial_bins(b(2)+1), &
          radial_bins(b(3)), radial_bins(b(3)+1), &
          radial_bins(b(4)), radial_bins(b(4)+1), &
          radial_bins(b(5)), radial_bins(b(5)+1), &
          radial_bins(b(6)), radial_bins(b(6)+1), &
          z_mean, z_sigma, o_mean, o_sigma, c_mean, c_sigma
      else
        write(unit_err, '(12(e14.7,1x),4(e14.7,1x))') &
          radial_bins(b(1)), radial_bins(b(1)+1), &
          radial_bins(b(2)), radial_bins(b(2)+1), &
          radial_bins(b(3)), radial_bins(b(3)+1), &
          radial_bins(b(4)), radial_bins(b(4)+1), &
          radial_bins(b(5)), radial_bins(b(5)+1), &
          radial_bins(b(6)), radial_bins(b(6)+1), &
          z_mean, z_sigma, c_mean, c_sigma
      end if
    end do
    close(unit_num)
    close(unit_err)
    deallocate(xi0_all)
    print '(" wrote ",i0," jackknife realisation rows to ",a)', &
      cfg%n_configs_4pcf * cfg%njk, trim(mode_output_file(trim(mode)))//'.jk'
  end subroutine write_4pcf_jackknife

  ! ---------------------------------------------------------------------------
  ! Cleanup direction lookup arrays
  ! ---------------------------------------------------------------------------
  subroutine cleanup_direction_lookup()
    if (allocated(dir_x)) deallocate(dir_x)
    if (allocated(px)) deallocate(px)
    if (allocated(py)) deallocate(py)
    if (allocated(pz)) deallocate(pz)
    if (allocated(dir_y)) deallocate(dir_y)
    if (allocated(dir_z)) deallocate(dir_z)
  end subroutine cleanup_direction_lookup

end module query_4pcf_module
