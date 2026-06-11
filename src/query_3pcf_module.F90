module query_3pcf_module
  use kdtree2_precision_module
  use iso_fortran_env, only: int8
  use config_module
  use graph_utils_module, only: find_dist, find_normal
  implicit none
contains

  ! Merge-walk implementation (default).
  ! For each hub i and each of its neighbors id1, performs a two-pointer walk
  ! through hub i's tail [k1+1..nn2] and id1's sorted neighbor list to find
  ! common neighbors in O(m) per k1 instead of O(m log m) binary searches.
  subroutine query_graph_3pcf_all(istart, iend)
    integer, intent(in) :: istart, iend
    integer :: i, k1, nn2, id1, nn_id1, a, b, bin
    integer(int8) :: ind1, ind2, ind3
    real(kdkind) :: wi_w1

    if (.not. cfg%half_graph) then
      print *, 'ERROR: merge-walk 3PCF requires half_graph=.true.'
      stop
    end if
    if (cfg%rank == 0) print *, 'Performing 3pcf (all configurations, merge-walk)'
    if (cfg%rank == 0) print *, 'begin querying the graph'

    !$OMP PARALLEL DO schedule(dynamic) &
    !$OMP& private(i, k1, nn2, id1, nn_id1, a, b, ind1, ind2, ind3, bin, wi_w1) &
    !$OMP& shared(weights, output, buffer, cfg, bintable) &
    !$OMP& reduction(+:N2) &
    !$OMP& reduction(+:N3)
    do i = istart, iend
      if (buffer(i) == 1) cycle

      nn2 = output(i)%nn

      do k1 = 1, nn2
        ind1   = output(i)%dist(k1)
        id1    = output(i)%id(k1)
        wi_w1  = weights(i) * weights(id1)
        nn_id1 = output(id1)%nn

        ! Two-pointer walk: hub's tail [k1+1..nn2] vs id1's full neighbor list.
        ! Both lists are sorted by node ID, so we advance the smaller pointer.
        ! Starting at k1+1 on the hub side preserves the k2>k1 uniqueness constraint
        ! since the list is ID-sorted (all hub tail entries have ID > id1).
        a = k1 + 1
        b = 1
        do while (a <= nn2 .and. b <= nn_id1)
          if (output(i)%id(a) == output(id1)%id(b)) then
            ind2 = output(i)%dist(a)
            ind3 = output(id1)%dist(b)
            if (ind2 /= 0 .and. ind3 /= 0) then
              bin = bintable(ind1, ind2, ind3, 1)
              if (cfg%RSD) then
                call find_normal(output(i)%mu(k1), output(i)%mu(a), ind2)
                if (i > cfg%num_data .and. id1 > cfg%num_data .and. &
                    output(i)%id(a) > cfg%num_data) then
                  N3(bin, ind2, 3) = N3(bin, ind2, 3) - wi_w1 * weights(output(i)%id(a))
                end if
                N2(bin, ind2, 3) = N2(bin, ind2, 3) + wi_w1 * weights(output(i)%id(a))
              else
                if (i > cfg%num_data .and. id1 > cfg%num_data .and. &
                    output(i)%id(a) > cfg%num_data) then
                  N3(bin, 1, 3) = N3(bin, 1, 3) - wi_w1 * weights(output(i)%id(a))
                end if
                N2(bin, 1, 3) = N2(bin, 1, 3) + wi_w1 * weights(output(i)%id(a))
              end if
            end if
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

    call write_3pcf_results()
  end subroutine query_graph_3pcf_all

  ! Original binary-search implementation, kept for benchmarking.
  ! Uses find_dist (O(log m) binary search) to check each candidate v2.
  subroutine query_graph_3pcf_all_bsearch(istart, iend)
    integer, intent(in) :: istart, iend
    integer :: i, k1, k2, nn2, id1, id2, bin
    integer(int8) :: ind1, ind2, ind3
    real(kdkind) :: wi_w1

    if (cfg%rank == 0) print *, 'Performing 3pcf (all configurations, binary-search)'
    if (cfg%rank == 0) print *, 'begin querying the graph'

    !$OMP PARALLEL DO schedule(dynamic) &
    !$OMP& private(i, k1, k2, nn2, id1, id2, ind1, ind2, ind3, bin, wi_w1) &
    !$OMP& shared(weights, output, buffer, cfg, bintable) &
    !$OMP& reduction(+:N2) &
    !$OMP& reduction(+:N3)
    do i = istart, iend
      if (buffer(i) == 1) cycle

      nn2 = output(i)%nn

      do k1 = 1, nn2
        ind1 = output(i)%dist(k1)
        id1 = output(i)%id(k1)
        wi_w1 = weights(i) * weights(id1)

        do k2 = k1 + 1, nn2
          ind2 = output(i)%dist(k2)
          if (ind2 == 0) cycle

          id2 = output(i)%id(k2)
          call find_dist(id1, id2, ind3)
          if (ind3 == 0) cycle

          bin = bintable(ind1, ind2, ind3, 1)

          if (cfg%RSD) then
            call find_normal(output(i)%mu(k1), output(i)%mu(k2), ind2)
            if (i > cfg%num_data .and. id1 > cfg%num_data .and. id2 > cfg%num_data) then
              N3(bin, ind2, 3) = N3(bin, ind2, 3) - wi_w1 * weights(id2)
            end if
            N2(bin, ind2, 3) = N2(bin, ind2, 3) + wi_w1 * weights(id2)
          else
            if (i > cfg%num_data .and. id1 > cfg%num_data .and. id2 > cfg%num_data) then
              N3(bin, 1, 3) = N3(bin, 1, 3) - wi_w1 * weights(id2)
            end if
            N2(bin, 1, 3) = N2(bin, 1, 3) + wi_w1 * weights(id2)
          end if
        end do
      end do
    end do
    !$OMP END PARALLEL DO

  end subroutine query_graph_3pcf_all_bsearch

  subroutine query_graph_equilateral_triangle(istart, iend)
    integer, intent(in) :: istart, iend
    integer :: i, k1, k2, nn2, id1, id2
    integer(int8) :: ind1, ind2, ind3

    if (cfg%rank == 0) print *, 'begin querying the graph'

    !$OMP PARALLEL DO schedule(dynamic) &
    !$OMP& private(i, k1, k2, nn2, id1, id2, ind1, ind2, ind3) &
    !$OMP& shared(weights, output, buffer, cfg) &
    !$OMP& reduction(+:N2, N3)
    do i = istart, iend
      if (buffer(i) == 1) cycle

      nn2 = output(i)%nn

      do k1 = 1, nn2
        ind1 = output(i)%dist(k1)
        id1 = output(i)%id(k1)

        do k2 = k1, nn2
          if (k2 == k1) cycle
          ind2 = output(i)%dist(k2)
          if (ind1 /= ind2) cycle

          id2 = output(i)%id(k2)
          call find_dist(id1, id2, ind3)
          if (ind3 /= ind1) cycle

          if (cfg%RSD) then
            call find_normal(output(i)%mu(k1), output(i)%mu(k2), ind2)
            N2(ind1, ind2, 3) = N2(ind1, ind2, 3) + weights(i) * weights(id1) * weights(id2)
          else
            if (i > cfg%num_data .and. id1 > cfg%num_data .and. id2 > cfg%num_data) then
              N3(ind1, 1, 3) = N3(ind1, 1, 3) - weights(i) * weights(id1) * weights(id2)
            end if
            N2(ind1, 1, 3) = N2(ind1, 1, 3) + weights(i) * weights(id1) * weights(id2)
          end if
        end do
      end do
    end do
    !$OMP END PARALLEL DO

    call write_equilateral_results()
  end subroutine query_graph_equilateral_triangle

  subroutine write_3pcf_results()
    integer :: i, j, k, l, unit_num

    if (cfg%rank /= cfg%master) return

    unit_num = 30
    open(unit_num, file=trim(cfg%output_file), status='unknown')
    if (cfg%RSD) then
      write(unit_num, *) 'r1 min, r1 max, r2 min, r2 max, r3 min, r3 max, mu min, mu max, NNN, RRR, 3pcf (zeta)'
    else
      write(unit_num, *) 'r1 min, r1 max, r2 min, r2 max, r3 min, r3 max, NNN, RRR, 3pcf (zeta)'
    end if

    do i = 1, cfg%nbins
      do j = i, cfg%nbins
        do k = j, cfg%nbins
          ! Read the pre-computed bin index (set by create_binlookup; no mutation here)
          associate(bin => bintable(i, j, k, 1))
          if (cfg%RSD) then
            do l = 1, cfg%nmu
              write(unit_num, '(11(e14.7,1x))') radial_bins(i), radial_bins(i+1), &
                radial_bins(j), radial_bins(j+1), radial_bins(k), radial_bins(k+1), &
                ((float(l)-1.)/cfg%mu_scale/2.), (float(l)/cfg%mu_scale/2.), &
                N2(bin, l, 3), N3(bin, l, 3), N2(bin, l, 3) / N3(bin, l, 3)
            end do
          else
            write(unit_num, '(9(e14.7,1x))') radial_bins(i), radial_bins(i+1), &
              radial_bins(j), radial_bins(j+1), radial_bins(k), radial_bins(k+1), &
              N2(bin, 1, 3), N3(bin, 1, 3), N2(bin, 1, 3) / N3(bin, 1, 3)
          end if
          end associate
        end do
      end do
    end do
    close(unit_num)
  end subroutine write_3pcf_results

  subroutine write_equilateral_results()
    integer :: l, k, unit_num

    if (cfg%rank /= cfg%master) return

    unit_num = 30
    open(unit_num, file=trim(cfg%output_file), status='unknown')
    do l = 1, cfg%nbins
      do k = 1, cfg%nmu
        write(unit_num, '(8(e14.7,1x))') radial_bins(l), radial_bins(l+1), &
             ((float(k)-1.)/cfg%mu_scale/2.), (float(k)/cfg%mu_scale/2.), &
             N2(l, k, 3), N3(l, k, 3), N2(l, k, 3) / N3(l, k, 3)
      end do
    end do
    close(unit_num)
  end subroutine write_equilateral_results

end module query_3pcf_module
