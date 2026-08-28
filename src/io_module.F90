module io_module
  use kdtree2_precision_module
  use config_module
  implicit none

  ! Catalogue lines: 'x y z weight' or 'x y z' (weight defaults to 1).
  ! Blank lines and lines starting with '#' are skipped; anything else that
  ! fails to parse is a hard error (a silent mis-parse would corrupt results).
  integer, parameter :: MAX_LINE_LEN = 4096

contains

  subroutine count_files_2()
    implicit none

    cfg%num_data = count_catalog(cfg%file1, 'data')
    cfg%num_rand = 0
    if (cfg%rancat) cfg%num_rand = count_catalog(cfg%file2, 'random')

    if (cfg%num_data == 0) then
      print *, 'ERROR: no data points found in ', trim(cfg%file1)
      stop 1
    end if
    if (cfg%rancat .and. cfg%num_rand == 0) then
      print *, 'ERROR: no random points found in ', trim(cfg%file2)
      stop 1
    end if

    if (cfg%rank == 0) print *, 'Preparing to read ', cfg%num_data, 'data points'
    if (cfg%rank == 0) print *, 'Preparing to read ', cfg%num_rand, 'random points'
  end subroutine count_files_2

  subroutine read_files_2()
    implicit none
    integer :: i

    call read_catalog(cfg%file1, 1, cfg%num_data, 'data')
    if (cfg%rancat) call read_catalog(cfg%file2, cfg%num_data + 1, cfg%num_rand, 'random')

    ! Normalize weights
    weights(1:cfg%num_data) = weights(1:cfg%num_data) / sum(weights(1:cfg%num_data))
    weights(cfg%num_data+1:cfg%num_data+cfg%num_rand) = &
      -1.0d0 * weights(cfg%num_data+1:cfg%num_data+cfg%num_rand) / &
      sum(weights(cfg%num_data+1:cfg%num_data+cfg%num_rand))

    ! Power sums of the normalized data weights, used by the analytic
    ! random counts (pair/triple/quadruple weight normalizations)
    cfg%sw2 = sum(weights(1:cfg%num_data)**2)
    cfg%sw3 = sum(weights(1:cfg%num_data)**3)
    cfg%sw4 = sum(weights(1:cfg%num_data)**4)

    ! Periodic box: wrap all coordinates into [0, L)
    if (cfg%periodic) then
      do i = 1, cfg%num_data + cfg%num_rand
        points(1, i) = modulo(points(1, i), cfg%boxsize(1))
        points(2, i) = modulo(points(2, i), cfg%boxsize(2))
        points(3, i) = modulo(points(3, i), cfg%boxsize(3))
      end do
    end if

    if (cfg%rank == 0) print *, 'Finished reading data file'
    if (cfg%rank == 0) print *, 'sum of weights: ', sum(weights)
  end subroutine read_files_2

  ! Number of data lines in a catalogue file (blanks/'#' comments skipped).
  function count_catalog(fname, label) result(nlines)
    character(len=*), intent(in) :: fname, label
    integer :: nlines
    integer :: ios, unit_num
    character(len=MAX_LINE_LEN) :: line

    open(newunit=unit_num, file=trim(fname), status='old', iostat=ios)
    if (ios /= 0) then
      print *, 'ERROR: cannot open ', label, ' file ', trim(fname)
      stop 1
    end if

    nlines = 0
    do
      read(unit_num, '(a)', iostat=ios) line
      if (ios /= 0) exit
      if (is_skippable(line)) cycle
      nlines = nlines + 1
    end do
    close(unit_num)
  end function count_catalog

  ! Read n catalogue points into points/weights/buffer starting at index i0.
  ! The column count (3 or 4) is detected from the first data line; with 3
  ! columns every weight is 1.  Reading is line-based, so a 3-column file can
  ! never consume two records per point (which would silently leave half the
  ! arrays uninitialized, as the old record-based reader did).
  subroutine read_catalog(fname, i0, n, label)
    character(len=*), intent(in) :: fname, label
    integer, intent(in) :: i0, n
    integer :: ios, unit_num, i, ncol, lineno
    character(len=MAX_LINE_LEN) :: line
    real(kdkind) :: v4(4)

    if (cfg%rank == 0) print *, 'opening ', trim(fname)
    open(newunit=unit_num, file=trim(fname), status='old', iostat=ios)
    if (ios /= 0) then
      print *, 'ERROR: cannot open ', label, ' file ', trim(fname)
      stop 1
    end if

    ncol = 0
    lineno = 0
    i = 0
    do while (i < n)
      read(unit_num, '(a)', iostat=ios) line
      lineno = lineno + 1
      if (ios /= 0) then
        print '("ERROR: unexpected end of ",a," file ",a," after ",i0," of ",i0," points")', &
          label, trim(fname), i, n
        stop 1
      end if
      if (is_skippable(line)) cycle

      if (ncol == 0) then
        ! Detect the column count from the first data line.
        read(line, *, iostat=ios) v4
        if (ios == 0) then
          ncol = 4
        else
          read(line, *, iostat=ios) v4(1:3)
          if (ios /= 0) then
            print '("ERROR: cannot parse line ",i0," of ",a," (expected x y z [weight]): ",a)', &
              lineno, trim(fname), trim(line)
            stop 1
          end if
          ncol = 3
          if (cfg%rank == 0) print *, 'No weight column in ', trim(fname), ': using weight = 1'
        end if
      else
        read(line, *, iostat=ios) v4(1:ncol)
        if (ios /= 0) then
          print '("ERROR: cannot parse line ",i0," of ",a,": ",a)', &
            lineno, trim(fname), trim(line)
          stop 1
        end if
      end if
      if (ncol == 3) v4(4) = 1.0d0

      i = i + 1
      points(1:3, i0 + i - 1) = v4(1:3)
      weights(i0 + i - 1) = v4(4)
    end do
    close(unit_num)
  end subroutine read_catalog

  pure function is_skippable(line) result(skip)
    character(len=*), intent(in) :: line
    logical :: skip
    character(len=len(line)) :: t
    t = adjustl(line)
    skip = (len_trim(t) == 0 .or. t(1:1) == '#')
  end function is_skippable

    ! Read delete-one jackknife region labels, 1..cfg%njk, one integer per line,
  ! in the same order as the corresponding catalogue. Points with a label
  ! outside 1..njk (or with no label file) are treated as belonging to no
  ! region, so they are never deleted -- that is the conservative choice.
  ! With -njk but no -jkgal/-jkran the labels are computed internally from
  ! the angular positions instead (assign_jk_regions_angular).
  subroutine read_jk_regions()
    integer :: i, ios, unit_num, r, nbad
    ! Allocated unconditionally: the OpenACC data clauses in the GPU kernel
    ! name `region` whether or not jackknife is requested, and mapping an
    ! unallocated allocatable is a runtime error.
    allocate(region(cfg%num_data + cfg%num_rand))
    region = 0
    if (cfg%njk <= 0) return
    if (cfg%jk_internal) then
      call assign_jk_regions_angular()
      return
    end if
    nbad = 0
    if (len_trim(cfg%jkgal) > 0) then
      open(newunit=unit_num, file=trim(cfg%jkgal), status='old', iostat=ios)
      if (ios /= 0) then
        print *, 'ERROR: cannot open -jkgal file ', trim(cfg%jkgal)
        stop
      end if
      do i = 1, cfg%num_data
        read(unit_num, *, iostat=ios) r
        if (ios /= 0) then
          print *, 'ERROR: -jkgal file has fewer rows than the galaxy catalogue'
          stop
        end if
        if (r < 1 .or. r > cfg%njk) then
          nbad = nbad + 1
        else
          region(i) = r
        end if
      end do
      close(unit_num)
    end if
    if (len_trim(cfg%jkran) > 0) then
      open(newunit=unit_num, file=trim(cfg%jkran), status='old', iostat=ios)
      if (ios /= 0) then
        print *, 'ERROR: cannot open -jkran file ', trim(cfg%jkran)
        stop
      end if
      do i = cfg%num_data + 1, cfg%num_data + cfg%num_rand
        read(unit_num, *, iostat=ios) r
        if (ios /= 0) then
          print *, 'ERROR: -jkran file has fewer rows than the random catalogue'
          stop
        end if
        if (r < 1 .or. r > cfg%njk) then
          nbad = nbad + 1
        else
          region(i) = r
        end if
      end do
      close(unit_num)
    end if
    if (cfg%rank == 0) then
      print *, 'read jackknife regions; unassigned points: ', &
               count(region == 0), ' out of range: ', nbad
    end if
  end subroutine read_jk_regions

  ! ---------------------------------------------------------------------------
  ! Internal jackknife partition on the SKY (angular coordinates only).
  !
  ! Regions are defined purely by the direction of each point as seen from
  ! the observer at the origin -- never by its radial distance -- so that a
  ! deleted region removes a full line-of-sight cone.  Angular and radial
  ! variations have different systematic sources in real data (masking,
  ! seeing, star contamination vs redshift errors, selection); a partition
  ! that mixed them would conflate the two in the error estimate.
  !
  ! Construction: nband ~ sqrt(njk) latitude bands split at equal-count
  ! quantiles of sin(dec) = z/r, each band cut into phi = atan2(y, x) slices
  ! at equal-count quantiles, with the slice counts summing to exactly njk.
  ! Boundaries are computed from the RANDOM catalogue (it traces the angular
  ! selection with the lowest shot noise) and applied identically to both
  ! catalogues, so data and randoms always share the same region footprint.
  ! The labels are written to <out>.jkgal / <out>.jkran in the exact format
  ! -jkgal/-jkran read back, for reproducibility.
  ! ---------------------------------------------------------------------------
  subroutine assign_jk_regions_angular()
    integer :: nsrc, i0src, nband, i, b, k, q, idx, ntot, unit_num
    integer, allocatable :: nslice(:), roff(:), cnt_d(:), cnt_r(:), nb_src(:)
    real(kdkind), allocatable :: s_src(:), phi_src(:), sedge(:), band_phi(:)
    real(kdkind), allocatable :: phedge(:,:)
    real(kdkind) :: r2, sv, phv
    real(kdkind), parameter :: PI = 3.141592653589793238462643383279502884197d0

    ntot = cfg%num_data + cfg%num_rand

    ! Boundary source: randoms (validated present by parseOptions).
    i0src = cfg%num_data + 1
    nsrc  = cfg%num_rand

    if (cfg%njk > nsrc) then
      print *, 'ERROR: -njk exceeds the number of random points'
      stop 1
    end if
    if (nsrc / cfg%njk < 100 .and. cfg%rank == 0) then
      print *, 'WARNING: fewer than 100 randoms per jackknife region; the'
      print *, '         region boundaries will be noisy -- reduce -njk'
    end if

    ! Angular coordinates of the boundary source.
    allocate(s_src(nsrc), phi_src(nsrc))
    do i = 1, nsrc
      call point_angles(points(:, i0src + i - 1), sv, phv)
      s_src(i) = sv
      phi_src(i) = phv
    end do

    ! Latitude bands: counts proportional to each band's slice count so every
    ! final region ends up with ~nsrc/njk source points.
    nband = max(1, nint(sqrt(real(cfg%njk, kdkind))))
    allocate(nslice(nband), roff(nband))
    nslice = cfg%njk / nband
    do b = 1, mod(cfg%njk, nband)
      nslice(b) = nslice(b) + 1
    end do
    roff(1) = 0
    do b = 2, nband
      roff(b) = roff(b - 1) + nslice(b - 1)
    end do

    ! Band edges at cumulative equal-count quantiles of sorted sin(dec).
    call hsort_real(s_src, nsrc)
    allocate(sedge(nband + 1))
    sedge(1) = -1.0001d0
    sedge(nband + 1) = 1.0001d0
    k = 0
    do b = 1, nband - 1
      k = k + nslice(b)
      idx = nint(real(nsrc, kdkind) * real(k, kdkind) / real(cfg%njk, kdkind))
      idx = min(max(idx, 1), nsrc - 1)
      sedge(b + 1) = 0.5d0 * (s_src(idx) + s_src(idx + 1))
    end do

    ! Per-band phi slice edges at equal-count quantiles of the band's phis.
    ! s_src is sorted, so recompute band membership from the original points.
    allocate(phedge(maxval(nslice) + 1, nband), band_phi(nsrc), nb_src(nband))
    nb_src = 0
    do b = 1, nband
      k = 0
      do i = 1, nsrc
        call point_angles(points(:, i0src + i - 1), sv, phv)
        if (sv > sedge(b) .and. sv <= sedge(b + 1)) then
          k = k + 1
          band_phi(k) = phv
        end if
      end do
      nb_src(b) = k
      if (k < nslice(b)) then
        print *, 'ERROR: jackknife latitude band ', b, ' holds ', k, &
                 ' randoms for ', nslice(b), ' slices; reduce -njk'
        stop 1
      end if
      call hsort_real(band_phi, k)
      phedge(1, b) = -PI - 0.001d0
      phedge(nslice(b) + 1, b) = PI + 0.001d0
      do q = 1, nslice(b) - 1
        idx = nint(real(k, kdkind) * real(q, kdkind) / real(nslice(b), kdkind))
        idx = min(max(idx, 1), k - 1)
        phedge(q + 1, b) = 0.5d0 * (band_phi(idx) + band_phi(idx + 1))
      end do
    end do

    ! Assign every point (data first, then randoms) with the SAME boundaries.
    allocate(cnt_d(cfg%njk), cnt_r(cfg%njk))
    cnt_d = 0
    cnt_r = 0
    do i = 1, ntot
      call point_angles(points(:, i), sv, phv)
      b = 1
      do while (b < nband .and. sv > sedge(b + 1))
        b = b + 1
      end do
      q = 1
      do while (q < nslice(b) .and. phv > phedge(q + 1, b))
        q = q + 1
      end do
      region(i) = roff(b) + q
      if (i <= cfg%num_data) then
        cnt_d(region(i)) = cnt_d(region(i)) + 1
      else
        cnt_r(region(i)) = cnt_r(region(i)) + 1
      end if
    end do

    if (cfg%rank == 0) then
      print '(" jackknife: ",i0," angular regions (",i0," sin(dec) bands)")', &
        cfg%njk, nband
      print '("   data per region:    min ",i0,"  max ",i0)', &
        minval(cnt_d), maxval(cnt_d)
      print '("   randoms per region: min ",i0,"  max ",i0)', &
        minval(cnt_r), maxval(cnt_r)
      if (any(cnt_d == 0)) &
        print *, 'WARNING: some jackknife regions contain no data points'
    end if

    ! Persist the labels in the -jkgal/-jkran input format, so the exact
    ! partition can be re-used (or inspected) with external label files.
    if (cfg%rank == cfg%master) then
      open(newunit=unit_num, file=trim(cfg%output_file)//'.jkgal', status='unknown')
      do i = 1, cfg%num_data
        write(unit_num, '(i0)') region(i)
      end do
      close(unit_num)
      open(newunit=unit_num, file=trim(cfg%output_file)//'.jkran', status='unknown')
      do i = cfg%num_data + 1, ntot
        write(unit_num, '(i0)') region(i)
      end do
      close(unit_num)
      print *, 'wrote jackknife region labels to ', &
        trim(cfg%output_file)//'.jkgal', ' and ', trim(cfg%output_file)//'.jkran'
    end if
  end subroutine assign_jk_regions_angular

  ! sin(dec) = z/r and phi = atan2(y, x) of one point as seen from the origin.
  ! A point exactly at the origin has no direction; park it at (0, 0).
  subroutine point_angles(p, s, ph)
    real(kdkind), intent(in) :: p(3)
    real(kdkind), intent(out) :: s, ph
    real(kdkind) :: r
    r = sqrt(p(1)*p(1) + p(2)*p(2) + p(3)*p(3))
    if (r > 0.0d0) then
      s = p(3) / r
      ph = atan2(p(2), p(1))
    else
      s = 0.0d0
      ph = 0.0d0
    end if
  end subroutine point_angles

  ! In-place ascending heapsort of a(1:n).  Self-contained so the partition
  ! has no dependency on the kd-tree internals.
  subroutine hsort_real(a, n)
    real(kdkind), intent(inout) :: a(:)
    integer, intent(in) :: n
    integer :: i, j, child
    real(kdkind) :: tmp
    do i = n / 2, 1, -1
      call sift_down(i, n)
    end do
    do i = n, 2, -1
      tmp = a(1)
      a(1) = a(i)
      a(i) = tmp
      call sift_down(1, i - 1)
    end do
  contains
    subroutine sift_down(start, last)
      integer, intent(in) :: start, last
      real(kdkind) :: v
      j = start
      v = a(j)
      do
        child = 2 * j
        if (child > last) exit
        if (child < last) then
          if (a(child + 1) > a(child)) child = child + 1
        end if
        if (a(child) <= v) exit
        a(j) = a(child)
        j = child
      end do
      a(j) = v
    end subroutine sift_down
  end subroutine hsort_real

end module io_module
