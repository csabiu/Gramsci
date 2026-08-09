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
    integer :: id3, jr1, jr2, jr3
    integer(int8) :: ind1, ind2, ind3
    real(kdkind) :: wi_w1, w3, wprod

    if (.not. cfg%half_graph) then
      print *, 'ERROR: merge-walk 3PCF requires half_graph=.true.'
      stop
    end if
    if (cfg%rank == 0) print *, 'Performing 3pcf (all configurations, merge-walk)'
    if (cfg%rank == 0) print *, 'begin querying the graph'

    !$OMP PARALLEL DO schedule(dynamic) &
    !$OMP& private(i, k1, nn2, id1, nn_id1, a, b, ind1, ind2, ind3, bin, wi_w1) &
    !$OMP& private(id3, jr1, jr2, jr3, w3, wprod) &
    !$OMP& shared(weights, output, buffer, cfg, bintable, region) &
    !$OMP& reduction(+:N2) &
    !$OMP& reduction(+:N3) &
    !$OMP& reduction(+:N2jk) &
    !$OMP& reduction(+:N3jk)
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
              id3 = output(i)%id(a)
              w3 = weights(id3)
              wprod = wi_w1 * w3
              ! Jackknife bookkeeping: accumulate this triplet against each
              ! DISTINCT region it touches. Deduplication matters -- a triplet
              ! with two nodes in the same region must be subtracted once, not
              ! twice, when that region is deleted.
              if (cfg%njk > 0) then
                jr1 = region(i)
                jr2 = region(id1)
                jr3 = region(id3)
              end if
              if (cfg%RSD) then
                call find_normal(output(i)%mu(k1), output(i)%mu(a), ind2)
                if (i > cfg%num_data .and. id1 > cfg%num_data .and. &
                    id3 > cfg%num_data) then
                  N3(bin, ind2, 3) = N3(bin, ind2, 3) - wprod
                  if (cfg%njk > 0) then
                    if (jr1 > 0) &
                      N3jk(bin, ind2, 3, jr1) = N3jk(bin, ind2, 3, jr1) - wprod
                    if (jr2 > 0 .and. jr2 /= jr1) &
                      N3jk(bin, ind2, 3, jr2) = N3jk(bin, ind2, 3, jr2) - wprod
                    if (jr3 > 0 .and. jr3 /= jr1 .and. jr3 /= jr2) &
                      N3jk(bin, ind2, 3, jr3) = N3jk(bin, ind2, 3, jr3) - wprod
                  end if
                end if
                N2(bin, ind2, 3) = N2(bin, ind2, 3) + wprod
                if (cfg%njk > 0) then
                  if (jr1 > 0) &
                    N2jk(bin, ind2, 3, jr1) = N2jk(bin, ind2, 3, jr1) + wprod
                  if (jr2 > 0 .and. jr2 /= jr1) &
                    N2jk(bin, ind2, 3, jr2) = N2jk(bin, ind2, 3, jr2) + wprod
                  if (jr3 > 0 .and. jr3 /= jr1 .and. jr3 /= jr2) &
                    N2jk(bin, ind2, 3, jr3) = N2jk(bin, ind2, 3, jr3) + wprod
                end if
              else
                if (i > cfg%num_data .and. id1 > cfg%num_data .and. &
                    id3 > cfg%num_data) then
                  N3(bin, 1, 3) = N3(bin, 1, 3) - wprod
                  if (cfg%njk > 0) then
                    if (jr1 > 0) &
                      N3jk(bin, 1, 3, jr1) = N3jk(bin, 1, 3, jr1) - wprod
                    if (jr2 > 0 .and. jr2 /= jr1) &
                      N3jk(bin, 1, 3, jr2) = N3jk(bin, 1, 3, jr2) - wprod
                    if (jr3 > 0 .and. jr3 /= jr1 .and. jr3 /= jr2) &
                      N3jk(bin, 1, 3, jr3) = N3jk(bin, 1, 3, jr3) - wprod
                  end if
                end if
                N2(bin, 1, 3) = N2(bin, 1, 3) + wprod
                if (cfg%njk > 0) then
                  if (jr1 > 0) &
                    N2jk(bin, 1, 3, jr1) = N2jk(bin, 1, 3, jr1) + wprod
                  if (jr2 > 0 .and. jr2 /= jr1) &
                    N2jk(bin, 1, 3, jr2) = N2jk(bin, 1, 3, jr2) + wprod
                  if (jr3 > 0 .and. jr3 /= jr1 .and. jr3 /= jr2) &
                    N2jk(bin, 1, 3, jr3) = N2jk(bin, 1, 3, jr3) + wprod
                end if
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
    call write_3pcf_jackknife()
  end subroutine query_graph_3pcf_all

  ! Original binary-search implementation, kept for benchmarking.
  ! Uses find_dist (O(log m) binary search) to check each candidate v2.
  subroutine query_graph_3pcf_all_bsearch(istart, iend)
    integer, intent(in) :: istart, iend
    integer :: i, k1, k2, nn2, id1, id2, bin
    integer(int8) :: ind1, ind2, ind3
    real(kdkind) :: wi_w1

    if (cfg%njk > 0) then
      print *, 'ERROR: -njk is implemented for the merge-walk 3PCF only, '// &
               'not the binary-search variant.'
      stop
    end if
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

    call write_3pcf_results()
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

  ! Write the delete-one jackknife realisations. Realisation m omits every
  ! triplet touching region m:  N_m = N_total - N_touching(m).
  ! Weight renormalisation is deliberately NOT applied: deleting a region
  ! scales the data and random weight sums by the same factor to the accuracy
  ! that the randoms trace the selection, and that factor cancels in the ratio
  ! zeta = N2/N3.
  subroutine write_3pcf_jackknife()
    integer :: i, j, k, l, m, unit_num
    real(kdkind) :: n2m, n3m

    if (cfg%rank /= cfg%master) return
    if (cfg%njk <= 0) return

    unit_num = 31
    open(unit_num, file=trim(cfg%output_file)//'.jk', status='unknown')
    write(unit_num, *) '# delete-one jackknife realisations, njk = ', cfg%njk
    if (cfg%RSD) then
      write(unit_num, *) '# r1min r1max r2min r2max r3min r3max mumin mumax ireal NNN RRR zeta'
    else
      write(unit_num, *) '# r1min r1max r2min r2max r3min r3max ireal NNN RRR zeta'
    end if

    do i = 1, cfg%nbins
      do j = i, cfg%nbins
        do k = j, cfg%nbins
          associate(bin => bintable(i, j, k, 1))
          do l = 1, cfg%nmu
            do m = 1, cfg%njk
              n2m = N2(bin, l, 3) - N2jk(bin, l, 3, m)
              n3m = N3(bin, l, 3) - N3jk(bin, l, 3, m)
              if (cfg%RSD) then
                write(unit_num, '(8(e14.7,1x),i6,1x,3(e14.7,1x))') &
                  radial_bins(i), radial_bins(i+1), radial_bins(j), radial_bins(j+1), &
                  radial_bins(k), radial_bins(k+1), &
                  ((float(l)-1.)/cfg%mu_scale/2.), (float(l)/cfg%mu_scale/2.), &
                  m, n2m, n3m, n2m / n3m
              else
                write(unit_num, '(6(e14.7,1x),i6,1x,3(e14.7,1x))') &
                  radial_bins(i), radial_bins(i+1), radial_bins(j), radial_bins(j+1), &
                  radial_bins(k), radial_bins(k+1), m, n2m, n3m, n2m / n3m
              end if
            end do
          end do
          end associate
        end do
      end do
    end do
    close(unit_num)
    print *, 'wrote jackknife realisations to ', trim(cfg%output_file)//'.jk'
  end subroutine write_3pcf_jackknife

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
