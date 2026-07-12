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

    unit_num = 20
    open(unit_num, file=trim(fname), status='old', iostat=ios)
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

    unit_num = 20
    if (cfg%rank == 0) print *, 'opening ', trim(fname)
    open(unit_num, file=trim(fname), status='old', iostat=ios)
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
      buffer(i0 + i - 1) = 0
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

end module io_module
