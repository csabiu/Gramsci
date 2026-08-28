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
    !$OMP& shared(weights, output, buffer, cfg, region) &
    !$OMP& reduction(+:N2, N3) &
    !$OMP& reduction(+:N2jk, N3jk)
    do i = istart, iend
      if (buffer(i) == 1) cycle

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

    unit_num = 30
    print *, 'writing output to: ', trim(mode_output_file('2pcf'))
    open(unit_num, file=trim(mode_output_file('2pcf')), status='unknown')
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
  end subroutine write_2pcf_results

  ! Write the delete-one jackknife realisations and their error summary.
  ! Realisation m omits every pair touching region m: N_m = N_total - N_jk(m).
  ! Weight renormalisation is deliberately NOT applied (see the 3PCF note:
  ! the deleted-region factor cancels in the ratio xi = N2/N3 to the accuracy
  ! that the randoms trace the selection).
  subroutine write_2pcf_jackknife()
    integer :: l, k, m, unit_num, unit_err
    real(kdkind) :: n2m, n3m, xi_m(cfg%njk), xi_mean, xi_sigma

    if (cfg%rank /= cfg%master) return
    if (cfg%njk <= 0) return

    unit_num = 31
    unit_err = 32
    open(unit_num, file=trim(mode_output_file('2pcf'))//'.jk', status='unknown')
    call write_provenance(unit_num)
    write(unit_num, '(a,i0)') '# delete-one jackknife realisations, njk = ', cfg%njk
    write(unit_num, '(a)') '# rmin rmax mumin mumax ireal NN RR xi'
    open(unit_err, file=trim(mode_output_file('2pcf'))//'.jkerr', status='unknown')
    call write_provenance(unit_err)
    write(unit_err, '(a,i0)') '# delete-one jackknife error, njk = ', cfg%njk
    write(unit_err, '(a)') '# rmin rmax mumin mumax xi_mean_jk xi_sigma_jk'

    do l = 1, cfg%nbins
      do k = 1, cfg%nmu
        do m = 1, cfg%njk
          n2m = N2(l, k, 3) - N2jk(l, k, 3, m)
          n3m = N3(l, k, 3) - N3jk(l, k, 3, m)
          if (n3m /= 0.0d0) then
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
      end do
    end do
    close(unit_num)
    close(unit_err)
    print *, 'wrote jackknife realisations to ', trim(mode_output_file('2pcf'))//'.jk'
  end subroutine write_2pcf_jackknife

end module query_2pcf_module
