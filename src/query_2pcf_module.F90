module query_2pcf_module
  use kdtree2_precision_module
  use iso_fortran_env, only: int8
  use config_module
  implicit none
contains

  subroutine query_graph_2pcf(istart, iend)
    integer, intent(in) :: istart, iend
    integer :: i, k1, nn2, id1
    integer :: jr1, jr2
    integer(int8) :: ind1, mu
    real(kdkind) :: wpair

    if (cfg%rank == 0) then
      print *, 'begin querying the graph'
      print *, 'number of mu bins:', cfg%nmu
    end if

    !$OMP PARALLEL DO schedule(dynamic) private(i, k1, nn2, id1, ind1, mu) &
    !$OMP& private(jr1, jr2, wpair) &
    !$OMP& shared(weights, output, cfg, region) &
    !$OMP& reduction(+:N2, N3) &
    !$OMP& reduction(+:N2jk, N3jk)
    do i = istart, iend
      nn2 = output(i)%nn

      do k1 = 1, nn2
        ind1 = output(i)%dist(k1)
        id1 = output(i)%id(k1)
        mu = 1
        if (cfg%RSD) mu = output(i)%mu(k1)
        wpair = weights(i) * weights(id1)

        ! Jackknife bookkeeping: accumulate this pair against each DISTINCT
        ! region it touches (a pair with both points in the same region must
        ! be subtracted once, not twice, when that region is deleted).
        if (cfg%njk > 0) then
          jr1 = region(i)
          jr2 = region(id1)
        end if

        if (i > cfg%num_data .and. id1 > cfg%num_data) then
          N3(ind1, mu, 3) = N3(ind1, mu, 3) + wpair
          if (cfg%njk > 0) then
            if (jr1 > 0) &
              N3jk(ind1, mu, 3, jr1) = N3jk(ind1, mu, 3, jr1) + wpair
            if (jr2 > 0 .and. jr2 /= jr1) &
              N3jk(ind1, mu, 3, jr2) = N3jk(ind1, mu, 3, jr2) + wpair
          end if
        end if
        N2(ind1, mu, 3) = N2(ind1, mu, 3) + wpair
        if (cfg%njk > 0) then
          if (jr1 > 0) &
            N2jk(ind1, mu, 3, jr1) = N2jk(ind1, mu, 3, jr1) + wpair
          if (jr2 > 0 .and. jr2 /= jr1) &
            N2jk(ind1, mu, 3, jr2) = N2jk(ind1, mu, 3, jr2) + wpair
        end if
      end do
    end do
    !$OMP END PARALLEL DO

    call write_2pcf_results()
    call write_2pcf_jackknife()
  end subroutine query_graph_2pcf

  subroutine write_2pcf_results()
    use analytic_randoms_module, only: analytic_rr
    integer :: l, k, unit_num
    real(kdkind) :: rr_an(cfg%nbins), xi

    if (cfg%rank /= cfg%master) return

    ! Analytic mode: RR from shell-volume fractions of the periodic box;
    ! NN holds pure DD, so the natural estimator xi = DD/RR - 1 applies
    ! (with a random catalogue the signed weights make NN/RR the LS
    ! estimator directly, no -1).
    if (cfg%analytic) then
      call analytic_rr(rr_an)
      do l = 1, cfg%nbins
        N3(l, 1:cfg%nmu, 3) = rr_an(l) / cfg%nmu
      end do
    end if

    print *, 'writing output to: ', trim(mode_output_file('2pcf'))
    open(newunit=unit_num, file=trim(mode_output_file('2pcf')), status='unknown')
    call write_provenance(unit_num)
    write(unit_num, '(a)') '# r min, r max, mu min, mu max, NN, RR, 2pcf (xi)'

    do l = 1, cfg%nbins
      do k = 1, cfg%nmu
        ! Empty bin (no RR pairs): write 0, not Inf/NaN — same convention
        ! as the 4PCF and jackknife writers.
        if (N3(l, k, 3) == 0.0d0) then
          xi = 0.0d0
        else if (cfg%analytic) then
          xi = N2(l, k, 3) / N3(l, k, 3) - 1.0d0
        else
          xi = N2(l, k, 3) / N3(l, k, 3)
        end if
        write(unit_num, '(8(e14.7,1x))') radial_bins(l), radial_bins(l+1), &
             ((float(k)-1.)/cfg%mu_scale)-1., (float(k)/cfg%mu_scale)-1.0, &
             N2(l, k, 3), N3(l, k, 3), xi
      end do
    end do
    close(unit_num)

    ! ---- Legendre multipoles (plane-parallel z in a box; midpoint LOS in
    ! survey mode with -nmu > 1): xi_ell = (2l+1) * S_ell / RR, with S_ell
    ! the pair Legendre sums accumulated at graph build.  Under the signed
    ! weights S_ell is already the correlation-only part for ell > 0 (the
    ! uniform term integrates to zero against L_ell); in analytic mode the
    ! monopole alone needs the "-1" (see compute_2pcf_for_4pcf).
    if (cfg%pair_mult) then
      block
        integer :: unit_m
        real(kdkind) :: rr_l, x0, x2, x4
        open(newunit=unit_m, file=trim(mode_output_file('2pcf'))//'.mult', &
             status='unknown')
        call write_provenance(unit_m)
        if (cfg%periodic) then
          write(unit_m, '(a)') '# 2PCF Legendre multipoles, plane-parallel line of sight (z axis)'
        else
          write(unit_m, '(a)') '# 2PCF Legendre multipoles, midpoint line of sight'
        end if
        write(unit_m, '(a)') '# r min, r max, xi0, xi2, xi4'
        do l = 1, cfg%nbins
          rr_l = sum(N3(l, 1:cfg%nmu, 3))
          if (rr_l /= 0.0d0) then
            if (cfg%analytic) then
              x0 = sum(N2(l, 1:cfg%nmu, 3)) / rr_l - 1.0d0
            else
              x0 = sum(N2(l, 1:cfg%nmu, 3)) / rr_l
            end if
            x2 = 5.0d0 * sum_pair_l2(l) / rr_l
            x4 = 9.0d0 * sum_pair_l4(l) / rr_l
          else
            x0 = 0.0d0
            x2 = 0.0d0
            x4 = 0.0d0
          end if
          write(unit_m, '(5(e14.7,1x))') radial_bins(l), radial_bins(l+1), &
            x0, x2, x4
        end do
        close(unit_m)
        print *, 'wrote 2PCF multipoles to ', trim(mode_output_file('2pcf'))//'.mult'
      end block
    end if
  end subroutine write_2pcf_results

  ! Write the delete-one jackknife realisations and their error summary.
  ! Realisation m omits every pair touching region m: N_m = N_total - N_jk(m).
  ! Weight renormalisation is deliberately NOT applied (see the 3PCF note:
  ! the deleted-region factor cancels in the ratio xi = N2/N3 to the accuracy
  ! that the randoms trace the selection).
  subroutine write_2pcf_jackknife()
    integer :: l, k, m, unit_num, unit_err, r
    real(kdkind) :: n2m, n3m, xi_m(cfg%njk), xi_mean, xi_sigma
    real(kdkind) :: xi_all(cfg%nbins * cfg%nmu, max(cfg%njk, 1))

    if (cfg%rank /= cfg%master) return
    if (cfg%njk <= 0) return

    open(newunit=unit_num, file=trim(mode_output_file('2pcf'))//'.jk', status='unknown')
    call write_provenance(unit_num)
    write(unit_num, '(a,i0)') '# delete-one jackknife realisations, njk = ', cfg%njk
    write(unit_num, '(a)') '# rmin rmax mumin mumax ireal NN RR xi'
    open(newunit=unit_err, file=trim(mode_output_file('2pcf'))//'.jkerr', status='unknown')
    call write_provenance(unit_err)
    write(unit_err, '(a,i0)') '# delete-one jackknife error, njk = ', cfg%njk
    write(unit_err, '(a)') '# rmin rmax mumin mumax xi_mean_jk xi_sigma_jk'

    do l = 1, cfg%nbins
      do k = 1, cfg%nmu
        do m = 1, cfg%njk
          n2m = N2(l, k, 3) - N2jk(l, k, 3, m)
          n3m = N3(l, k, 3) - N3jk(l, k, 3, m)
          ! Relative emptiness test, not an exact zero (JK_DENOM_TOL note)
          if (abs(n3m) > JK_DENOM_TOL * abs(N3(l, k, 3))) then
            xi_m(m) = n2m / n3m
          else
            xi_m(m) = 0.0d0
          end if
          write(unit_num, '(4(e14.7,1x),i6,1x,3(e14.7,1x))') &
            radial_bins(l), radial_bins(l+1), &
            ((float(k)-1.)/cfg%mu_scale)-1., (float(k)/cfg%mu_scale)-1.0, &
            m, n2m, n3m, xi_m(m)
        end do
        call jk_mean_sigma(xi_m, xi_mean, xi_sigma)
        write(unit_err, '(6(e14.7,1x))') radial_bins(l), radial_bins(l+1), &
          ((float(k)-1.)/cfg%mu_scale)-1., (float(k)/cfg%mu_scale)-1.0, &
          xi_mean, xi_sigma
        r = (l - 1) * cfg%nmu + k
        xi_all(r, :) = xi_m
      end do
    end do
    close(unit_num)
    close(unit_err)
    call write_jk_covariance(trim(mode_output_file('2pcf'))//'.jkcov', xi_all)
    print *, 'wrote jackknife realisations to ', trim(mode_output_file('2pcf'))//'.jk'
  end subroutine write_2pcf_jackknife

end module query_2pcf_module
