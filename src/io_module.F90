module io_module
  use kdtree2_precision_module
  use config_module
  implicit none
contains

  subroutine count_files()
    implicit none
    integer :: ios, unit_num
    real(kdkind) :: aux

    cfg%num_data = 0
    if (cfg%DOMPI) cfg%file1 = trim(str(cfg%rank+1)) // '.loadnodes'

    unit_num = 20
    print *, trim(cfg%file1)
    open(unit_num, file=trim(cfg%file1), status='old', iostat=ios)
    if (ios /= 0) then
      print *, 'ERROR: cannot open file ', trim(cfg%file1)
      stop
    end if

    do
      read(unit_num, *, iostat=ios) aux
      if (ios /= 0) exit
      cfg%num_data = cfg%num_data + 1
    end do
    close(unit_num)

    if (cfg%rank == 0) print *, 'Preparing to read ', cfg%num_data, 'data points'
  end subroutine count_files

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

  subroutine read_files()
    implicit none
    integer :: i, ios, unit_num

    unit_num = 20
    if (cfg%rank == 0) print *, 'opening ', trim(cfg%file1)
    open(unit_num, file=trim(cfg%file1), status='old', iostat=ios)
    if (ios /= 0) then
      print *, 'ERROR: cannot open file ', trim(cfg%file1)
      stop
    end if

    do i = 1, cfg%num_data
      if (cfg%cut) then
        read(unit_num, *, iostat=ios) points(1:3, i), weights(i), buffer(i)
      else
        read(unit_num, *, iostat=ios) points(1:3, i), weights(i)
        buffer(i) = 0
      end if
      if (ios /= 0) exit
    end do
    close(unit_num)

    if (cfg%rank == 0) print *, 'Finished reading data file'
    if (cfg%rank == 0) print *, 'there are ', cfg%num_data - sum(buffer), ' data points inside buffer'
  end subroutine read_files

  ! Read delete-one jackknife region labels, 1..cfg%njk, one integer per line,
  ! in the same order as the corresponding catalogue. Points with a label
  ! outside 1..njk (or with no label file) are treated as belonging to no
  ! region, so they are never deleted -- that is the conservative choice.
  subroutine read_jk_regions()
    integer :: i, ios, unit_num, r, nbad
    ! Allocated unconditionally: the OpenACC data clauses in the GPU kernel
    ! name `region` whether or not jackknife is requested, and mapping an
    ! unallocated allocatable is a runtime error.
    allocate(region(cfg%num_data + cfg%num_rand))
    region = 0
    if (cfg%njk <= 0) return
    nbad = 0
    unit_num = 41
    if (len_trim(cfg%jkgal) > 0) then
      open(unit_num, file=trim(cfg%jkgal), status='old', iostat=ios)
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
      open(unit_num, file=trim(cfg%jkran), status='old', iostat=ios)
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

    if (cfg%rank == 0) print *, 'Finished reading data file'
    if (cfg%rank == 0) print *, 'sum of weights: ', sum(weights)
  end subroutine read_files_2

end module io_module
