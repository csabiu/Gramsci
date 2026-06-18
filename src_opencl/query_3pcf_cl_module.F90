! ---------------------------------------------------------------------------
! OpenCL-offloaded 3PCF queries (all configurations + equilateral).
!
! Mirrors src_gpu/query_3pcf_gpu_module.F90 but dispatches to the OpenCL
! kernels k_3pcf_all / k_3pcf_equi.  Single-pass only: the whole CSR graph plus
! accumulators are uploaded once.  Each work-item owns one accumulator column
! (no atomics); the ngang columns are summed into N2/N3 in double on the host.
!
! Limitation: RSD (nmu > 1) is not offloaded — the caller routes those to the
! CPU implementation, exactly as the OpenACC build does.
! ---------------------------------------------------------------------------
module query_3pcf_cl_module
  use iso_c_binding
  use iso_fortran_env, only: int8, int32, int64, real32
  use kdtree2_precision_module
  use config_module
  use csr_cl_module
  use cl_env_module
  use query_3pcf_module, only: write_3pcf_results, write_equilateral_results
  implicit none
  private
  public :: query_graph_3pcf_cl, query_graph_equilateral_cl

contains

  ! Choose the number of work-items / accumulator columns ("gangs").
  ! The Apple GPU needs many in-flight work-items to hide memory latency, so
  ! keep this large (throughput); the interleaved-bucket launcher seeds its
  ! first window from the hubs/work-item ratio to stay occupied.  Fewer columns
  ! also keeps each fp32 partial small (good accuracy after the double host
  ! reduction).  Bounded by the device's max single-buffer size.
  function pick_ngang(cb, n_hubs) result(ngang)
    integer, intent(in) :: cb
    integer(int64), intent(in) :: n_hubs
    integer(int64) :: ngang, cap_cols
    ngang = min(n_hubs, 65536_int64)
    ! keep each partial buffer (cb*ngang*4 bytes) within max alloc and <= 1 GB
    cap_cols = min(cl_max_alloc, 1073741824_int64) / (int(cb,int64) * 4_int64)
    ngang = min(ngang, cap_cols)
    ngang = max(ngang, 1_int64)
  end function pick_ngang

  ! Pack the float weights once (config_module stores them as double).
  subroutine pack_weights(wf)
    real(real32), allocatable, intent(out) :: wf(:)
    integer(int64) :: n, i
    n = cfg%num_data + cfg%num_rand
    allocate(wf(n))
    do i = 1, n
      wf(i) = real(weights(i), real32)
    end do
  end subroutine pack_weights

  subroutine query_graph_3pcf_cl(istart, iend)
    integer, intent(in) :: istart, iend
    integer :: cb, nb, bin, g
    integer(int64) :: ngang, n_hubs, ncol
    real(real32), allocatable :: wf(:), hn(:), hr(:)
    integer(c_intptr_t) :: b_ptr, b_id, b_dist, b_w, b_buf, b_bt3, b_pn, b_pr, kern
    real(kdkind) :: acc

    if (cfg%RSD) then
      print *, 'ERROR: query_graph_3pcf_cl called with RSD mode; route to CPU'
      stop
    end if
    if (cfg%rank == 0) print *, 'Performing 3PCF (all configurations, OpenCL bsearch)'

    nb = cfg%nbins
    cb = cfg%config_bins
    n_hubs = cfg%num_data + cfg%num_rand
    ngang = pick_ngang(cb, n_hubs)
    ncol = int(cb,int64) * ngang
    if (cfg%rank == 0) print '("  ngang=",i0,"  config_bins=",i0)', ngang, cb

    call pack_weights(wf)
    allocate(hn(ncol), hr(ncol))
    hn = 0.0_real32
    hr = 0.0_real32

    ! ---- Upload ----
    b_ptr  = cl_buf_in_i64(csr_ptr,  int(size(csr_ptr),int64))
    b_id   = cl_buf_in_i32(csr_id,   csr_total_edges)
    b_dist = cl_buf_in_i8 (csr_dist, csr_total_edges)
    b_w    = cl_buf_in_f32(wf,       n_hubs)
    b_buf  = cl_buf_in_i32(buffer,   n_hubs)
    ! bintable(:,:,:,1) is stored column-major == the kernel's flatten index,
    ! so its first nbins^3 elements ARE the lookup table -> upload directly.
    b_bt3  = cl_buf_in_i32(bintable,  int(nb,int64)**3)
    b_pn   = cl_buf_zeroed_f32(hn,   ncol)
    b_pr   = cl_buf_zeroed_f32(hr,   ncol)

    ! ---- Launch ----
    kern = cl_kernel_get('k_3pcf_all')
    call cl_arg_mem(kern, 0, b_ptr)
    call cl_arg_mem(kern, 1, b_id)
    call cl_arg_mem(kern, 2, b_dist)
    call cl_arg_mem(kern, 3, b_w)
    call cl_arg_mem(kern, 4, b_buf)
    call cl_arg_mem(kern, 5, b_bt3)
    call cl_arg_mem(kern, 6, b_pn)
    call cl_arg_mem(kern, 7, b_pr)
    call cl_arg_i32(kern, 8,  nb)
    call cl_arg_i32(kern, 9,  cb)
    call cl_arg_i32(kern, 10, cfg%num_data)
    call cl_arg_i32(kern, 11, istart)
    call cl_arg_i32(kern, 12, iend)
    call cl_arg_i32(kern, 13, int(ngang, int32))
    if (.not. cl_run_complete(kern, 14, 15, 16, 17, ngang)) then
      if (cfg%rank == 0) print *, '  single GPU launch hit the watchdog; retrying tiled (slower)'
      call cl_write_f32(b_pn, hn, ncol)   ! hn/hr are still all zeros here
      call cl_write_f32(b_pr, hr, ncol)
      call cl_run_bucketed(kern, 14, 15, 16, 17, ngang, n_hubs)
    end if

    ! ---- Read back + reduce columns in double ----
    call cl_read_f32(b_pn, hn, ncol)
    call cl_read_f32(b_pr, hr, ncol)
    do bin = 1, cb
      acc = 0.0d0
      do g = 0, int(ngang) - 1
        acc = acc + real(hn(int(g,int64)*cb + bin), kdkind)
      end do
      N2(bin, 1, 3) = N2(bin, 1, 3) + acc
      acc = 0.0d0
      do g = 0, int(ngang) - 1
        acc = acc + real(hr(int(g,int64)*cb + bin), kdkind)
      end do
      N3(bin, 1, 3) = N3(bin, 1, 3) + acc
    end do

    call cl_release(b_ptr);  call cl_release(b_id);  call cl_release(b_dist)
    call cl_release(b_w);    call cl_release(b_buf); call cl_release(b_bt3)
    call cl_release(b_pn);   call cl_release(b_pr);  call cl_release_kernel(kern)
    deallocate(wf, hn, hr)

    call write_3pcf_results()
  end subroutine query_graph_3pcf_cl

  subroutine query_graph_equilateral_cl(istart, iend)
    integer, intent(in) :: istart, iend
    integer :: nb, bin, g
    integer(int64) :: ngang, n_hubs, ncol
    real(real32), allocatable :: wf(:), hn(:), hr(:)
    integer(c_intptr_t) :: b_ptr, b_id, b_dist, b_w, b_buf, b_pn, b_pr, kern
    real(kdkind) :: acc

    if (cfg%RSD) then
      print *, 'ERROR: query_graph_equilateral_cl called with RSD mode; route to CPU'
      stop
    end if
    if (cfg%rank == 0) print *, 'Performing equilateral 3PCF (OpenCL bsearch)'

    nb = cfg%nbins
    n_hubs = cfg%num_data + cfg%num_rand
    ngang = pick_ngang(nb, n_hubs)
    ncol = int(nb,int64) * ngang
    if (cfg%rank == 0) print '("  ngang=",i0,"  nbins=",i0)', ngang, nb

    call pack_weights(wf)
    allocate(hn(ncol), hr(ncol))
    hn = 0.0_real32
    hr = 0.0_real32

    b_ptr  = cl_buf_in_i64(csr_ptr,  int(size(csr_ptr),int64))
    b_id   = cl_buf_in_i32(csr_id,   csr_total_edges)
    b_dist = cl_buf_in_i8 (csr_dist, csr_total_edges)
    b_w    = cl_buf_in_f32(wf,       n_hubs)
    b_buf  = cl_buf_in_i32(buffer,   n_hubs)
    b_pn   = cl_buf_zeroed_f32(hn,   ncol)
    b_pr   = cl_buf_zeroed_f32(hr,   ncol)

    kern = cl_kernel_get('k_3pcf_equi')
    call cl_arg_mem(kern, 0, b_ptr)
    call cl_arg_mem(kern, 1, b_id)
    call cl_arg_mem(kern, 2, b_dist)
    call cl_arg_mem(kern, 3, b_w)
    call cl_arg_mem(kern, 4, b_buf)
    call cl_arg_mem(kern, 5, b_pn)
    call cl_arg_mem(kern, 6, b_pr)
    call cl_arg_i32(kern, 7,  nb)
    call cl_arg_i32(kern, 8,  nb)
    call cl_arg_i32(kern, 9,  cfg%num_data)
    call cl_arg_i32(kern, 10, istart)
    call cl_arg_i32(kern, 11, iend)
    call cl_arg_i32(kern, 12, int(ngang, int32))
    if (.not. cl_run_complete(kern, 13, 14, 15, 16, ngang)) then
      if (cfg%rank == 0) print *, '  single GPU launch hit the watchdog; retrying tiled (slower)'
      call cl_write_f32(b_pn, hn, ncol)
      call cl_write_f32(b_pr, hr, ncol)
      call cl_run_bucketed(kern, 13, 14, 15, 16, ngang, n_hubs)
    end if

    call cl_read_f32(b_pn, hn, ncol)
    call cl_read_f32(b_pr, hr, ncol)
    do bin = 1, nb
      acc = 0.0d0
      do g = 0, int(ngang) - 1
        acc = acc + real(hn(int(g,int64)*nb + bin), kdkind)
      end do
      N2(bin, 1, 3) = N2(bin, 1, 3) + acc
      acc = 0.0d0
      do g = 0, int(ngang) - 1
        acc = acc + real(hr(int(g,int64)*nb + bin), kdkind)
      end do
      N3(bin, 1, 3) = N3(bin, 1, 3) + acc
    end do

    call cl_release(b_ptr);  call cl_release(b_id);  call cl_release(b_dist)
    call cl_release(b_w);    call cl_release(b_buf)
    call cl_release(b_pn);   call cl_release(b_pr);  call cl_release_kernel(kern)
    deallocate(wf, hn, hr)

    call write_equilateral_results()
  end subroutine query_graph_equilateral_cl

end module query_3pcf_cl_module
