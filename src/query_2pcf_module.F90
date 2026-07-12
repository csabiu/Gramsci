module query_2pcf_module
  use kdtree2_precision_module
  use iso_fortran_env, only: int8
  use config_module
  implicit none
contains

  subroutine query_graph_2pcf(istart, iend)
    integer, intent(in) :: istart, iend
    integer :: i, k1, nn2, id1
    integer(int8) :: ind1, mu

    if (cfg%rank == 0) then
      print *, 'begin querying the graph'
      print *, 'number of mu bins:', cfg%nmu
    end if

    !$OMP PARALLEL DO schedule(dynamic) private(i, k1, nn2, id1, ind1, mu) &
    !$OMP& shared(weights, output, buffer, cfg) &
    !$OMP& reduction(+:N2, N3)
    do i = istart, iend
      if (buffer(i) == 1) cycle

      nn2 = output(i)%nn

      do k1 = 1, nn2
        ind1 = output(i)%dist(k1)
        id1 = output(i)%id(k1)
        mu = 1
        if (cfg%RSD) mu = output(i)%mu(k1)

        if (i > cfg%num_data .and. id1 > cfg%num_data) then
          N3(ind1, mu, 3) = N3(ind1, mu, 3) + weights(i) * weights(id1)
        end if
        N2(ind1, mu, 3) = N2(ind1, mu, 3) + weights(i) * weights(id1)
      end do
    end do
    !$OMP END PARALLEL DO

    call write_2pcf_results()
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
    print *, 'writing output to: ', trim(cfg%output_file)
    open(unit_num, file=trim(cfg%output_file), status='unknown')
    write(unit_num, *) 'r min, r max, mu min, mu max, NN, RR, 2pcf (xi)'

    do l = 1, cfg%nbins
      do k = 1, cfg%nmu
        if (cfg%analytic) then
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

end module query_2pcf_module
