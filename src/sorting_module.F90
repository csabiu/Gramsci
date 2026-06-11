module sorting_module
  use iso_fortran_env, only: int8
  implicit none
contains

  ! Insertion sort for three parallel arrays (id, dist, mu) keyed on id.
  ! O(n) for nearly-sorted data (common with KD-tree neighbor results),
  ! low overhead for small arrays typical of neighbor counts (10-200).
  subroutine sort2(x, y, z, n)
    integer, intent(inout) :: x(:)
    integer(int8), intent(inout) :: y(:), z(:)
    integer, intent(in) :: n
    integer :: i, j, key_x
    integer(int8) :: key_y, key_z

    do i = 2, n
      key_x = x(i)
      key_y = y(i)
      key_z = z(i)
      j = i - 1
      do while (j >= 1 .and. x(j) > key_x)
        x(j+1) = x(j)
        y(j+1) = y(j)
        z(j+1) = z(j)
        j = j - 1
      end do
      x(j+1) = key_x
      y(j+1) = key_y
      z(j+1) = key_z
    end do
  end subroutine sort2

  ! Insertion sort for two parallel arrays (id, dist) keyed on id.
  subroutine sort3(x, y, n)
    integer, intent(inout) :: x(:)
    integer(int8), intent(inout) :: y(:)
    integer, intent(in) :: n
    integer :: i, j, key_x
    integer(int8) :: key_y

    do i = 2, n
      key_x = x(i)
      key_y = y(i)
      j = i - 1
      do while (j >= 1 .and. x(j) > key_x)
        x(j+1) = x(j)
        y(j+1) = y(j)
        j = j - 1
      end do
      x(j+1) = key_x
      y(j+1) = key_y
    end do
  end subroutine sort3

  ! Insertion sort for four parallel arrays (id, dist, mu, phi) keyed on id.
  ! Used when direction pixel is stored for 4PCF parity computation.
  subroutine sort2_with_phi(x, y, z, p, n)
    integer, intent(inout) :: x(:)
    integer(int8), intent(inout) :: y(:), z(:), p(:)
    integer, intent(in) :: n
    integer :: i, j, key_x
    integer(int8) :: key_y, key_z, key_p

    do i = 2, n
      key_x = x(i)
      key_y = y(i)
      key_z = z(i)
      key_p = p(i)
      j = i - 1
      do while (j >= 1 .and. x(j) > key_x)
        x(j+1) = x(j)
        y(j+1) = y(j)
        z(j+1) = z(j)
        p(j+1) = p(j)
        j = j - 1
      end do
      x(j+1) = key_x
      y(j+1) = key_y
      z(j+1) = key_z
      p(j+1) = key_p
    end do
  end subroutine sort2_with_phi

end module sorting_module
