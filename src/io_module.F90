module io_module
  use kdtree2_precision_module
  use config_module
  implicit none
contains

  subroutine count_files_2()
    implicit none
    integer :: ios, unit_num
    real(kdkind) :: aux

    cfg%num_data = 0
    cfg%num_rand = 0
    unit_num = 20

    open(unit_num, file=trim(cfg%file1), status='old', iostat=ios)
    if (ios /= 0) then
      print *, 'ERROR: cannot open data file ', trim(cfg%file1)
      stop
    end if
    do
      read(unit_num, *, iostat=ios) aux
      if (ios /= 0) exit
      cfg%num_data = cfg%num_data + 1
    end do
    close(unit_num)

    if (cfg%rancat) then
      open(unit_num, file=trim(cfg%file2), status='old', iostat=ios)
      if (ios /= 0) then
        print *, 'ERROR: cannot open random file ', trim(cfg%file2)
        stop
      end if
      do
        read(unit_num, *, iostat=ios) aux
        if (ios /= 0) exit
        cfg%num_rand = cfg%num_rand + 1
      end do
      close(unit_num)
    end if

    if (cfg%rank == 0) print *, 'Preparing to read ', cfg%num_data, 'data points'
    if (cfg%rank == 0) print *, 'Preparing to read ', cfg%num_rand, 'random points'
  end subroutine count_files_2

  subroutine read_files_2()
    implicit none
    integer :: i, ios, unit_num

    unit_num = 20
    if (cfg%rank == 0) print *, 'opening ', trim(cfg%file1)
    open(unit_num, file=trim(cfg%file1), status='old', iostat=ios)
    if (ios /= 0) then
      print *, 'ERROR: cannot open data file ', trim(cfg%file1)
      stop
    end if

    do i = 1, cfg%num_data
      read(unit_num, *, iostat=ios) points(1:3, i), weights(i)
      if (ios /= 0) then
        print *, 'WARNING: read error at line ', i
        exit
      end if
      buffer(i) = 0
    end do
    close(unit_num)

    if (cfg%rancat) then
      if (cfg%rank == 0) print *, 'opening ', trim(cfg%file2)
      open(unit_num, file=trim(cfg%file2), status='old', iostat=ios)
      if (ios /= 0) then
        print *, 'ERROR: cannot open random file ', trim(cfg%file2)
        stop
      end if

      do i = cfg%num_data + 1, cfg%num_data + cfg%num_rand
        read(unit_num, *, iostat=ios) points(1:3, i), weights(i)
        if (ios /= 0) then
          print *, 'WARNING: read error at line ', i
          exit
        end if
        buffer(i) = 0
      end do
      close(unit_num)
    end if

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

end module io_module
