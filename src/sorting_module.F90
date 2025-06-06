module sorting_module
  use iso_fortran_env, only: int8
  implicit none
contains
  integer function findminimum(x, start, end)
    integer, intent(in) :: x(:)
    integer, intent(in) :: start, end
    integer :: minimum, location, i
    minimum = x(start)
    location = start
    do i = start + 1, end
       if (x(i) < minimum) then
          minimum = x(i)
          location = i
       end if
    end do
    findminimum = location
  end function findminimum

  subroutine swap(a, b)
    integer, intent(inout) :: a, b
    integer :: temp
    temp = a
    a = b
    b = temp
  end subroutine swap

  subroutine swap2(a, b)
    integer(int8), intent(inout) :: a, b
    integer(int8) :: temp
    temp = a
    a = b
    b = temp
  end subroutine swap2

  subroutine sort2(x, y, z, size)
    integer, intent(inout) :: x(:)
    integer(int8), intent(inout) :: y(:), z(:)
    integer, intent(in) :: size
    integer :: i, location
    do i = 1, size - 1
       location = findminimum(x, i, size)
       call swap(x(i), x(location))
       call swap2(y(i), y(location))
       call swap2(z(i), z(location))
    end do
  end subroutine sort2

  subroutine sort3(x, y, size)
    integer, intent(inout) :: x(:)
    integer(int8), intent(inout) :: y(:)
    integer, intent(in) :: size
    integer :: i, location
    do i = 1, size - 1
       location = findminimum(x, i, size)
       call swap(x(i), x(location))
       call swap2(y(i), y(location))
    end do
  end subroutine sort3
end module sorting_module
