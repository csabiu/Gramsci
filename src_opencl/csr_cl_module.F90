! ---------------------------------------------------------------------------
! CSR (Compressed Sparse Row) adjacency representation for the OpenCL backend.
!
! Same flat layout as the OpenACC build's csr_module, but portable: it does NOT
! use glibc malloc_trim (absent on macOS) or the OpenACC runtime.  The jagged
! graph output(:) is flattened into contiguous arrays that map directly onto
! OpenCL device buffers:
!
!   csr_ptr(i) .. csr_ptr(i+1)-1   — range of edge indices for node i (1-based)
!   csr_id(e)                       — neighbor node index for edge e
!   csr_dist(e)                     — distance-bin index (int8)
!   csr_phi(e)                      — direction-pixel index (int16, 4PCFp only:
!                                     up to ntheta*nphi = 32767 pixels, 256 on
!                                     the default 8x32 grid — too big for int8)
!
! Each jagged row is freed as it is copied so peak host RAM stays near one
! graph copy.  csr_ptr / csr_total_edges are 64-bit (graphs may exceed 2^31
! edges); node ids stay int32.
!
!   call build_csr()           after create_graph()
!   call deallocate_csr()      at cleanup
! ---------------------------------------------------------------------------
module csr_cl_module
  use iso_fortran_env, only: int8, int16, int32, int64
  use config_module
  implicit none

  integer(int64), allocatable :: csr_ptr(:)  ! shape (N+1), 1-based offsets
  integer(int32), allocatable :: csr_id(:)   ! shape (total_edges)
  integer(int8),  allocatable :: csr_dist(:) ! shape (total_edges)
  integer(int16), allocatable :: csr_phi(:)  ! shape (total_edges), 4PCFp only
  integer(int64) :: csr_total_edges = 0

contains

  subroutine build_csr()
    integer :: i, k, N
    integer(int64) :: base

    N = cfg%num_data + cfg%num_rand

    ! Pass 1: prefix sums -> csr_ptr (64-bit; can exceed 2^31-1 edges)
    allocate(csr_ptr(N + 1))
    csr_ptr(1) = 1
    do i = 1, N
      csr_ptr(i + 1) = csr_ptr(i) + output(i)%nn
    end do
    csr_total_edges = csr_ptr(N + 1) - 1

    allocate(csr_id(csr_total_edges))
    allocate(csr_dist(csr_total_edges))
    if (cfg%four_pcf_parity) allocate(csr_phi(csr_total_edges))

    ! Pass 2: copy from jagged output(:), freeing each row as we go
    do i = 1, N
      base = csr_ptr(i) - 1
      do k = 1, output(i)%nn
        csr_id  (base + k) = output(i)%id  (k)
        csr_dist(base + k) = output(i)%dist(k)
      end do
      if (cfg%four_pcf_parity .and. allocated(output(i)%phi)) then
        do k = 1, output(i)%nn
          csr_phi(base + k) = output(i)%phi(k)
        end do
      end if
      call output(i)%destroy()
    end do
    deallocate(output)

    if (cfg%rank == 0) &
      print '("CSR graph: ",i0," nodes, ",i0," edges")', N, csr_total_edges
  end subroutine build_csr

  ! Longest adjacency row in the CSR graph (sizes the per-work-item lmat scratch).
  function csr_max_row_len() result(mx)
    integer :: mx, i, n
    mx = 0
    n = cfg%num_data + cfg%num_rand
    do i = 1, n
      mx = max(mx, int(csr_ptr(i + 1) - csr_ptr(i)))
    end do
  end function csr_max_row_len

  subroutine deallocate_csr()
    if (allocated(csr_ptr))  deallocate(csr_ptr)
    if (allocated(csr_id))   deallocate(csr_id)
    if (allocated(csr_dist)) deallocate(csr_dist)
    if (allocated(csr_phi))  deallocate(csr_phi)
    csr_total_edges = 0
  end subroutine deallocate_csr

end module csr_cl_module
