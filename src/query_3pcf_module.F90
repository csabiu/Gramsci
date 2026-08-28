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
      stop 1
    end if
    if (cfg%rank == 0) print *, 'Performing 3pcf (all configurations, merge-walk)'
    if (cfg%rank == 0) print *, 'begin querying the graph'

    !$OMP PARALLEL DO schedule(dynamic) &
    !$OMP& private(i, k1, nn2, id1, nn_id1, a, b, ind1, ind2, ind3, bin, wi_w1) &
    !$OMP& private(id3, jr1, jr2, jr3, w3, wprod) &
    !$OMP& shared(weights, output, cfg, bintable, region) &
    !$OMP& reduction(+:N2) &
    !$OMP& reduction(+:N3) &
    !$OMP& reduction(+:N2jk) &
    !$OMP& reduction(+:N3jk)
    do i = istart, iend
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
    !$OMP& shared(weights, output, cfg, bintable) &
    !$OMP& reduction(+:N2) &
    !$OMP& reduction(+:N3)
    do i = istart, iend
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
    integer :: jr1, jr2, jr3
    integer(int8) :: ind1, ind2, ind3
    real(kdkind) :: wprod

    if (cfg%rank == 0) print *, 'begin querying the graph'

    !$OMP PARALLEL DO schedule(dynamic) &
    !$OMP& private(i, k1, k2, nn2, id1, id2, ind1, ind2, ind3) &
    !$OMP& private(jr1, jr2, jr3, wprod) &
    !$OMP& shared(weights, output, cfg, region) &
    !$OMP& reduction(+:N2, N3) &
    !$OMP& reduction(+:N2jk, N3jk)
    do i = istart, iend
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

          wprod = weights(i) * weights(id1) * weights(id2)
          ! Jackknife: accumulate against each DISTINCT region touched
          ! (see the all-configurations 3PCF for the dedup rationale).
          if (cfg%njk > 0) then
            jr1 = region(i)
            jr2 = region(id1)
            jr3 = region(id2)
          end if

          if (cfg%RSD) then
            call find_normal(output(i)%mu(k1), output(i)%mu(k2), ind2)
            if (i > cfg%num_data .and. id1 > cfg%num_data .and. id2 > cfg%num_data) then
              N3(ind1, ind2, 3) = N3(ind1, ind2, 3) - wprod
              if (cfg%njk > 0) then
                if (jr1 > 0) &
                  N3jk(ind1, ind2, 3, jr1) = N3jk(ind1, ind2, 3, jr1) - wprod
                if (jr2 > 0 .and. jr2 /= jr1) &
                  N3jk(ind1, ind2, 3, jr2) = N3jk(ind1, ind2, 3, jr2) - wprod
                if (jr3 > 0 .and. jr3 /= jr1 .and. jr3 /= jr2) &
                  N3jk(ind1, ind2, 3, jr3) = N3jk(ind1, ind2, 3, jr3) - wprod
              end if
            end if
            N2(ind1, ind2, 3) = N2(ind1, ind2, 3) + wprod
            if (cfg%njk > 0) then
              if (jr1 > 0) &
                N2jk(ind1, ind2, 3, jr1) = N2jk(ind1, ind2, 3, jr1) + wprod
              if (jr2 > 0 .and. jr2 /= jr1) &
                N2jk(ind1, ind2, 3, jr2) = N2jk(ind1, ind2, 3, jr2) + wprod
              if (jr3 > 0 .and. jr3 /= jr1 .and. jr3 /= jr2) &
                N2jk(ind1, ind2, 3, jr3) = N2jk(ind1, ind2, 3, jr3) + wprod
            end if
          else
            if (i > cfg%num_data .and. id1 > cfg%num_data .and. id2 > cfg%num_data) then
              N3(ind1, 1, 3) = N3(ind1, 1, 3) - wprod
              if (cfg%njk > 0) then
                if (jr1 > 0) &
                  N3jk(ind1, 1, 3, jr1) = N3jk(ind1, 1, 3, jr1) - wprod
                if (jr2 > 0 .and. jr2 /= jr1) &
                  N3jk(ind1, 1, 3, jr2) = N3jk(ind1, 1, 3, jr2) - wprod
                if (jr3 > 0 .and. jr3 /= jr1 .and. jr3 /= jr2) &
                  N3jk(ind1, 1, 3, jr3) = N3jk(ind1, 1, 3, jr3) - wprod
              end if
            end if
            N2(ind1, 1, 3) = N2(ind1, 1, 3) + wprod
            if (cfg%njk > 0) then
              if (jr1 > 0) &
                N2jk(ind1, 1, 3, jr1) = N2jk(ind1, 1, 3, jr1) + wprod
              if (jr2 > 0 .and. jr2 /= jr1) &
                N2jk(ind1, 1, 3, jr2) = N2jk(ind1, 1, 3, jr2) + wprod
              if (jr3 > 0 .and. jr3 /= jr1 .and. jr3 /= jr2) &
                N2jk(ind1, 1, 3, jr3) = N2jk(ind1, 1, 3, jr3) + wprod
            end if
          end if
        end do
      end do
    end do
    !$OMP END PARALLEL DO

    call write_equilateral_results()
    call write_equilateral_jackknife()
  end subroutine query_graph_equilateral_triangle

  subroutine write_3pcf_results()
    use analytic_randoms_module, only: analytic_rrr
    integer :: i, j, k, l, unit_num
    real(kdkind) :: zeta

    if (cfg%rank /= cfg%master) return

    ! Analytic mode: RRR from the periodic-box triangle kernel; NNN holds
    ! pure DDD, so DDD/RRR = 1 + xi(r1) + xi(r2) + xi(r3) + zeta and the
    ! connected 3PCF requires subtracting the internal 2PCF of the sides
    ! (with a random catalogue the signed weights give zeta = NNN/RRR
    ! directly).  xi_2pcf is filled by compute_2pcf_for_4pcf beforehand.
    if (cfg%analytic) then
      call analytic_rrr(N3(:, 1, 3))
    end if

    ! Internal pass (analytic 4PCF): store the connected 3PCF per config for
    ! the 4PCF subtraction instead of writing an output file.
    if (cfg%internal_3pcf) then
      if (.not. allocated(zeta3_internal)) allocate(zeta3_internal(cfg%config_bins))
      zeta3_internal = 0.0d0
      do i = 1, cfg%nbins
        do j = i, cfg%nbins
          do k = j, cfg%nbins
            associate(bin => bintable(i, j, k, 1))
            if (N3(bin, 1, 3) > 0.0d0) then
              zeta3_internal(bin) = N2(bin, 1, 3) / N3(bin, 1, 3) - 1.0d0 &
                                    - xi_2pcf(i) - xi_2pcf(j) - xi_2pcf(k)
            end if
            end associate
          end do
        end do
      end do
      return
    end if

    open(newunit=unit_num, file=trim(mode_output_file('3pcf')), status='unknown')
    call write_provenance(unit_num)
    if (cfg%RSD) then
      write(unit_num, '(a)') '# r1 min, r1 max, r2 min, r2 max, r3 min, r3 max, mu min, mu max, NNN, RRR, 3pcf (zeta)'
    else
      write(unit_num, '(a)') '# r1 min, r1 max, r2 min, r2 max, r3 min, r3 max, NNN, RRR, 3pcf (zeta)'
    end if

    do i = 1, cfg%nbins
      do j = i, cfg%nbins
        do k = j, cfg%nbins
          ! Read the pre-computed bin index (set by create_binlookup; no mutation here)
          associate(bin => bintable(i, j, k, 1))
          if (cfg%RSD) then
            do l = 1, cfg%nmu
              ! Empty bin (no RRR triplets): write 0, not Inf/NaN.
              if (N3(bin, l, 3) /= 0.0d0) then
                zeta = N2(bin, l, 3) / N3(bin, l, 3)
              else
                zeta = 0.0d0
              end if
              write(unit_num, '(11(e14.7,1x))') radial_bins(i), radial_bins(i+1), &
                radial_bins(j), radial_bins(j+1), radial_bins(k), radial_bins(k+1), &
                ((float(l)-1.)/cfg%mu_scale/2.), (float(l)/cfg%mu_scale/2.), &
                N2(bin, l, 3), N3(bin, l, 3), zeta
            end do
          else
            if (cfg%analytic) then
              if (N3(bin, 1, 3) > 0.0d0) then
                zeta = N2(bin, 1, 3) / N3(bin, 1, 3) - 1.0d0 &
                       - xi_2pcf(i) - xi_2pcf(j) - xi_2pcf(k)
              else
                zeta = 0.0d0
              end if
            else if (N3(bin, 1, 3) /= 0.0d0) then
              zeta = N2(bin, 1, 3) / N3(bin, 1, 3)
            else
              zeta = 0.0d0
            end if
            write(unit_num, '(9(e14.7,1x))') radial_bins(i), radial_bins(i+1), &
              radial_bins(j), radial_bins(j+1), radial_bins(k), radial_bins(k+1), &
              N2(bin, 1, 3), N3(bin, 1, 3), zeta
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
    integer :: i, j, k, l, m, unit_num, unit_err, r
    real(kdkind) :: n2m, n3m, zeta_m(cfg%njk), z_mean, z_sigma
    real(kdkind) :: z_all(cfg%config_bins * cfg%nmu, max(cfg%njk, 1))

    if (cfg%rank /= cfg%master) return
    if (cfg%njk <= 0) return
    r = 0
    ! The analytic 4PCF's internal 3PCF pass must never emit jackknife files
    ! (analytic mode rejects -njk anyway; this is belt and braces).
    if (cfg%internal_3pcf) return

    open(newunit=unit_num, file=trim(mode_output_file('3pcf'))//'.jk', status='unknown')
    call write_provenance(unit_num)
    write(unit_num, '(a,i0)') '# delete-one jackknife realisations, njk = ', cfg%njk
    if (cfg%RSD) then
      write(unit_num, '(a)') '# r1min r1max r2min r2max r3min r3max mumin mumax ireal NNN RRR zeta'
    else
      write(unit_num, '(a)') '# r1min r1max r2min r2max r3min r3max ireal NNN RRR zeta'
    end if
    open(newunit=unit_err, file=trim(mode_output_file('3pcf'))//'.jkerr', status='unknown')
    call write_provenance(unit_err)
    write(unit_err, '(a,i0)') '# delete-one jackknife error, njk = ', cfg%njk
    if (cfg%RSD) then
      write(unit_err, '(a)') '# r1min r1max r2min r2max r3min r3max mumin mumax zeta_mean_jk zeta_sigma_jk'
    else
      write(unit_err, '(a)') '# r1min r1max r2min r2max r3min r3max zeta_mean_jk zeta_sigma_jk'
    end if

    do i = 1, cfg%nbins
      do j = i, cfg%nbins
        do k = j, cfg%nbins
          associate(bin => bintable(i, j, k, 1))
          do l = 1, cfg%nmu
            do m = 1, cfg%njk
              n2m = N2(bin, l, 3) - N2jk(bin, l, 3, m)
              n3m = N3(bin, l, 3) - N3jk(bin, l, 3, m)
              if (n3m /= 0.0d0) then
                zeta_m(m) = n2m / n3m
              else
                zeta_m(m) = 0.0d0
              end if
              if (cfg%RSD) then
                write(unit_num, '(8(e14.7,1x),i6,1x,3(e14.7,1x))') &
                  radial_bins(i), radial_bins(i+1), radial_bins(j), radial_bins(j+1), &
                  radial_bins(k), radial_bins(k+1), &
                  ((float(l)-1.)/cfg%mu_scale/2.), (float(l)/cfg%mu_scale/2.), &
                  m, n2m, n3m, zeta_m(m)
              else
                write(unit_num, '(6(e14.7,1x),i6,1x,3(e14.7,1x))') &
                  radial_bins(i), radial_bins(i+1), radial_bins(j), radial_bins(j+1), &
                  radial_bins(k), radial_bins(k+1), m, n2m, n3m, zeta_m(m)
              end if
            end do
            call jk_mean_sigma(zeta_m, z_mean, z_sigma)
            if (cfg%RSD) then
              write(unit_err, '(10(e14.7,1x))') &
                radial_bins(i), radial_bins(i+1), radial_bins(j), radial_bins(j+1), &
                radial_bins(k), radial_bins(k+1), &
                ((float(l)-1.)/cfg%mu_scale/2.), (float(l)/cfg%mu_scale/2.), &
                z_mean, z_sigma
            else
              write(unit_err, '(8(e14.7,1x))') &
                radial_bins(i), radial_bins(i+1), radial_bins(j), radial_bins(j+1), &
                radial_bins(k), radial_bins(k+1), z_mean, z_sigma
            end if
            r = r + 1
            z_all(r, :) = zeta_m
          end do
          end associate
        end do
      end do
    end do
    close(unit_num)
    close(unit_err)
    call write_jk_covariance(trim(mode_output_file('3pcf'))//'.jkcov', z_all)
    print *, 'wrote jackknife realisations to ', trim(mode_output_file('3pcf'))//'.jk'
  end subroutine write_3pcf_jackknife

  subroutine write_equilateral_results()
    use analytic_randoms_module, only: analytic_rrr_equilateral
    integer :: l, k, unit_num
    real(kdkind) :: zeta

    if (cfg%rank /= cfg%master) return

    ! Analytic mode: equilateral configs (l,l,l); see write_3pcf_results.
    if (cfg%analytic) then
      call analytic_rrr_equilateral(N3(:, 1, 3))
    end if

    open(newunit=unit_num, file=trim(mode_output_file('equi')), status='unknown')
    call write_provenance(unit_num)
    write(unit_num, '(a)') '# r min, r max, mu min, mu max, NNN, RRR, equilateral 3pcf (zeta)'
    do l = 1, cfg%nbins
      do k = 1, cfg%nmu
        if (cfg%analytic) then
          if (N3(l, k, 3) > 0.0d0) then
            zeta = N2(l, k, 3) / N3(l, k, 3) - 1.0d0 - 3.0d0 * xi_2pcf(l)
          else
            zeta = 0.0d0
          end if
        else if (N3(l, k, 3) /= 0.0d0) then
          zeta = N2(l, k, 3) / N3(l, k, 3)
        else
          ! Empty bin (no RRR triplets): write 0, not Inf/NaN.
          zeta = 0.0d0
        end if
        write(unit_num, '(8(e14.7,1x))') radial_bins(l), radial_bins(l+1), &
             ((float(k)-1.)/cfg%mu_scale/2.), (float(k)/cfg%mu_scale/2.), &
             N2(l, k, 3), N3(l, k, 3), zeta
      end do
    end do
    close(unit_num)
  end subroutine write_equilateral_results

  ! Delete-one jackknife for the equilateral 3PCF; same conventions as
  ! write_3pcf_jackknife (no weight renormalisation, N_m = N - N_touching(m)).
  subroutine write_equilateral_jackknife()
    integer :: l, k, m, unit_num, unit_err, r
    real(kdkind) :: n2m, n3m, zeta_m(cfg%njk), z_mean, z_sigma
    real(kdkind) :: z_all(cfg%nbins * cfg%nmu, max(cfg%njk, 1))

    if (cfg%rank /= cfg%master) return
    if (cfg%njk <= 0) return

    open(newunit=unit_num, file=trim(mode_output_file('equi'))//'.jk', status='unknown')
    call write_provenance(unit_num)
    write(unit_num, '(a,i0)') '# delete-one jackknife realisations, njk = ', cfg%njk
    write(unit_num, '(a)') '# rmin rmax mumin mumax ireal NNN RRR zeta'
    open(newunit=unit_err, file=trim(mode_output_file('equi'))//'.jkerr', status='unknown')
    call write_provenance(unit_err)
    write(unit_err, '(a,i0)') '# delete-one jackknife error, njk = ', cfg%njk
    write(unit_err, '(a)') '# rmin rmax mumin mumax zeta_mean_jk zeta_sigma_jk'

    do l = 1, cfg%nbins
      do k = 1, cfg%nmu
        do m = 1, cfg%njk
          n2m = N2(l, k, 3) - N2jk(l, k, 3, m)
          n3m = N3(l, k, 3) - N3jk(l, k, 3, m)
          if (n3m /= 0.0d0) then
            zeta_m(m) = n2m / n3m
          else
            zeta_m(m) = 0.0d0
          end if
          write(unit_num, '(4(e14.7,1x),i6,1x,3(e14.7,1x))') &
            radial_bins(l), radial_bins(l+1), &
            ((float(k)-1.)/cfg%mu_scale/2.), (float(k)/cfg%mu_scale/2.), &
            m, n2m, n3m, zeta_m(m)
        end do
        call jk_mean_sigma(zeta_m, z_mean, z_sigma)
        write(unit_err, '(6(e14.7,1x))') radial_bins(l), radial_bins(l+1), &
          ((float(k)-1.)/cfg%mu_scale/2.), (float(k)/cfg%mu_scale/2.), &
          z_mean, z_sigma
        r = (l - 1) * cfg%nmu + k
        z_all(r, :) = zeta_m
      end do
    end do
    close(unit_num)
    close(unit_err)
    call write_jk_covariance(trim(mode_output_file('equi'))//'.jkcov', z_all)
    print *, 'wrote jackknife realisations to ', trim(mode_output_file('equi'))//'.jk'
  end subroutine write_equilateral_jackknife

end module query_3pcf_module
