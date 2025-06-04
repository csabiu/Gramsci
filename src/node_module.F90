module node_module
  use iso_fortran_env, only: int32, int8
  implicit none
  type :: node
    integer(kind=int32) :: nn = 0
    integer(kind=int32), allocatable :: id(:)
    integer(kind=int8), allocatable :: dist(:), mu(:)
  contains
    procedure :: init => node_init
    procedure :: destroy => node_destroy
  end type node
contains
  subroutine node_init(self, nsize)
    class(node), intent(inout) :: self
    integer, intent(in) :: nsize
    self%nn = nsize
    allocate(self%id(nsize))
    allocate(self%dist(nsize))
    allocate(self%mu(nsize))
  end subroutine node_init

  subroutine node_destroy(self)
    class(node), intent(inout) :: self
    if (allocated(self%id)) deallocate(self%id)
    if (allocated(self%dist)) deallocate(self%dist)
    if (allocated(self%mu)) deallocate(self%mu)
    self%nn = 0
  end subroutine node_destroy
end module node_module
