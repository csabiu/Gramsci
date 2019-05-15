module extension
  !
  ! defines in F90 some C commands
  ! These extensions are not completely standard in all F90/95 compilers
  !
  ! getEnvironment   : emulates getenv
  ! getArgument      : emulates getarg
  ! nArguments       : emulates iargc
!  IMPLICIT none

  integer(kind=4), private :: arg_shift = 0
!VF  integer(kind=I4B), private :: arg_shift = 1

  private
  public :: getArgument, nArguments, exit_with_status

  contains

#if (defined (GFORTRAN))

    ! ===========================================================
    function iargc ()
    ! ===========================================================
       integer iargc
       ! ===========================================================

       iargc=command_argument_count()
    end function

    ! ===========================================================
    subroutine getarg(num, res)
    ! ===========================================================
       integer, intent(in) :: num
       character(len=*), intent(out) :: res
       integer l, err
       ! ===========================================================
       call get_command_argument(num,res,l,err)
    end subroutine

#endif

    ! ===========================================================
    function nArguments() result(narg)
      ! ===========================================================
      integer(kind=4) :: narg
      ! ===========================================================

      narg = iargc() - arg_shift
!VF      narg = nargs() - arg_shift

      return
    end function nArguments
    ! ===========================================================
    ! ===========================================================
    subroutine getArgument(index, argument)
      ! ===========================================================
      integer(kind=4), intent(in) :: index
      character(len=*), intent(out) :: argument
      integer(kind=4) :: i1
      ! ===========================================================
      i1 = index + arg_shift
      call getarg(i1, argument)

      return
    end subroutine getArgument

    ! ===========================================================
    subroutine exit_with_status (code, msg)
      ! ===========================================================
      integer, intent(in) :: code
      character (len=*), intent(in), optional :: msg
      ! ===========================================================

      if (present(msg)) print *,trim(msg)
      print *,'program exits with exit code ', code

      call exit (code)

    end subroutine exit_with_status


end module extension
