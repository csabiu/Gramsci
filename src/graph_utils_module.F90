module graph_utils_module
  use kdtree2_precision_module
  use iso_fortran_env, only: int8
  use config_module, only: output, cfg
  implicit none
contains

  subroutine find_dist(id1, id2, ind)
    ! Binary search for node id2 in the sorted neighbor list of node id1.
    ! Returns the distance bin index in ind, or 0 if not found.
    implicit none
    integer, intent(in) :: id1, id2
    integer(int8), intent(out) :: ind
    integer :: lo, hi, mid

    if (id1 <= 0) then
      print *, 'ERROR: invalid node index', id1
      ind = 0
      return
    end if

    ! Quick bounds check
    if (output(id1)%nn == 0) then
      ind = 0
      return
    end if
    if (output(id1)%id(1) > id2 .or. output(id1)%id(output(id1)%nn) < id2) then
      ind = 0
      return
    end if

    ! Binary search
    lo = 1
    hi = output(id1)%nn
    do while (lo <= hi)
      mid = lo + (hi - lo) / 2
      if (output(id1)%id(mid) == id2) then
        ind = output(id1)%dist(mid)
        return
      else if (output(id1)%id(mid) < id2) then
        lo = mid + 1
      else
        hi = mid - 1
      end if
    end do

    ind = 0  ! not found
  end subroutine find_dist

  subroutine find_dist2(id1, id2, ind)
    ! Linear search fallback for node id2 in the neighbor list of node id1.
    implicit none
    integer, intent(in) :: id1, id2
    integer(int8), intent(out) :: ind
    integer :: i

    if (output(id1)%nn == 0) then
      ind = 0
      return
    end if
    if (output(id1)%id(1) > id2 .or. output(id1)%id(output(id1)%nn) < id2) then
      ind = 0
      return
    end if

    ind = 0
    do i = 1, output(id1)%nn
      if (output(id1)%id(i) == id2) then
        ind = output(id1)%dist(i)
        return
      end if
    end do
  end subroutine find_dist2

  subroutine find_normal(mu1, mu2, mun)
    implicit none
    integer(int8), intent(in) :: mu1, mu2
    integer(int8), intent(out) :: mun
    real :: mu11, mu22

    mu11 = ((mu1 - 0.5) / cfg%mu_scale) - 1.0
    mu22 = ((mu2 - 0.5) / cfg%mu_scale) - 1.0

    mun = int(floor((1.1547 * (0.75 - mu11*mu11 - mu22*mu22 + (mu11*mu22))**0.5) * cfg%nmu) + 1, int8)
    if (mun < 1) mun = 1
    if (mun > cfg%nmu) mun = int(cfg%nmu, int8)
  end subroutine find_normal

  pure function euclidean_dist(v1, v2) result(d)
    real(kdkind), intent(in) :: v1(3), v2(3)
    real(kdkind) :: d
    d = sqrt((v1(1)-v2(1))**2 + (v1(2)-v2(2))**2 + (v1(3)-v2(3))**2)
  end function euclidean_dist

end module graph_utils_module
