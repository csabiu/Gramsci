module node_module
  use iso_fortran_env, only: int32, int8, int16
  implicit none
  type :: node
    integer(kind=int32) :: nn = 0
    integer(kind=int32), allocatable :: id(:)
    integer(kind=int8), allocatable :: dist(:), mu(:)
    ! Direction pixel index (combined theta/phi bin), allocated only for 4PCF
    ! parity.  int16 so grids finer than 127 pixels are possible.
    integer(kind=int16), allocatable :: phi(:)
  contains
    procedure :: init => node_init
    procedure :: init_with_phi => node_init_with_phi
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

  subroutine node_init_with_phi(self, nsize)
    class(node), intent(inout) :: self
    integer, intent(in) :: nsize
    self%nn = nsize
    allocate(self%id(nsize))
    allocate(self%dist(nsize))
    allocate(self%mu(nsize))
    allocate(self%phi(nsize))
  end subroutine node_init_with_phi

  subroutine node_destroy(self)
    class(node), intent(inout) :: self
    if (allocated(self%id)) deallocate(self%id)
    if (allocated(self%dist)) deallocate(self%dist)
    if (allocated(self%mu)) deallocate(self%mu)
    if (allocated(self%phi)) deallocate(self%phi)
    self%nn = 0
  end subroutine node_destroy
end module node_module
