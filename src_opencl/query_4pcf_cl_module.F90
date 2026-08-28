! ---------------------------------------------------------------------------
! OpenCL-offloaded 4PCF queries (all configurations + parity decomposition).
!
! Pure binary-search strategy (matches the CPU query_graph_4pcf_bsearch): for
! each k1<k2<k3 the three closing edges are found with a device binary search.
!
! Unlike the OpenACC build, this does NOT precompute a per-hub connectivity
! matrix (lmat); lmat needs a per-work-item READ-WRITE global scratch buffer,
! whereas the bsearch kernel uses only read-only inputs plus a private
! accumulator column (structurally identical to the reliable 3PCF kernel) and
! is deterministic.  Each work-item owns one accumulator column (no atomics);
! the ngang columns are summed into N4/R4 in double on the host.
!
! Launch strategy (cl_run_complete + cl_run_bucketed): a single NDRange over
! all hubs is by far the fastest on Apple's GPU, but a long one trips the GPU
! watchdog, after which Apple's OpenCL returns PARTIAL results with no error and
! silently drops counts.  So we run ONE launch first and check per-work-item
! completion flags; only if it was truncated do we re-zero and re-run as several
! shorter interleaved-bucket launches (slower, but correct).
!
! Parity sign: the chirality sign_V is precomputed on the host (in double, with
! the same VOL_DEGEN_TOL as the CPU code) into a signv[ndir^3] table and looked
! up in the kernel — so the parity decomposition is bit-exact with the CPU
! reference even though the kernels run in fp32.
!
! Single-pass only: requires the CSR graph, the nbins^6 config table, and the
! accumulators to fit in device memory (max single buffer = cl_max_alloc).
! ---------------------------------------------------------------------------
module query_4pcf_cl_module
  use iso_c_binding
  use iso_fortran_env, only: int8, int32, int64, real32, real64
  use kdtree2_precision_module
  use config_module
  use csr_cl_module
  use cl_env_module
  use query_4pcf_module, only: write_4pcf_results, write_4pcf_results_noparity, &
                               dir_x, dir_y, dir_z, VOL_DEGEN_TOL
  implicit none
  private
  public :: query_graph_4pcf_cl, query_graph_4pcf_parity_cl

contains

  ! Size ngang (work-items / accumulator columns) from the device.  Each of the
  ! n_part accumulator buffers is ncfg*ngang floats (n_part counts the Kahan
  ! compensation buffers too); honour the single-buffer limit and keep the
  ! total within ~70% of device memory.
  function pick_ngang_4pcf(ncfg, n_part) result(ngang)
    integer, intent(in) :: ncfg, n_part
    integer(int64) :: ngang, part_cap, n_hubs, total
    n_hubs = cfg%num_data + cfg%num_rand
    part_cap = min(cl_max_alloc, 536870912_int64)
    ngang = max(1_int64, part_cap / (int(ncfg, int64) * 4_int64))
    ngang = min(ngang, n_hubs, 16384_int64)
    do
      total = int(n_part, int64) * ngang * int(ncfg, int64) * 4_int64
      if (total <= (cl_global_mem * 7_int64) / 10_int64 .or. ngang <= 64_int64) exit
      ngang = ngang / 2_int64
    end do
    ngang = max(ngang, 1_int64)
  end function pick_ngang_4pcf

  subroutine pack_weights(wf)
    real(real32), allocatable, intent(out) :: wf(:)
    integer(int64) :: n, i
    n = cfg%num_data + cfg%num_rand
    allocate(wf(n))
    do i = 1, n
      wf(i) = real(weights(i), real32)
    end do
  end subroutine pack_weights

  ! Guard: the nbins^6 config table must fit in one device buffer.
  subroutine check_bt6_fits(nb)
    integer, intent(in) :: nb
    integer(int64) :: bytes
    bytes = int(nb, int64)**6 * 4_int64
    if (bytes > cl_max_alloc) then
      print '("ERROR: nbins^6 config table is ",f6.2," GB > device max buffer ",f6.2," GB.")', &
            real(bytes)/1.0e9, real(cl_max_alloc)/1.0e9
      print *, '       Reduce -nbins, or use the NVIDIA OpenACC build (src_gpu).'
      stop 1
    end if
  end subroutine check_bt6_fits

  ! -------------------------------------------------------------------------
  ! 4PCF, all configurations (no parity).
  ! -------------------------------------------------------------------------
  subroutine query_graph_4pcf_cl(istart, iend)
    integer, intent(in) :: istart, iend
    integer :: nb, ncfg, c, g
    integer(int64) :: ngang, n_hubs, ncol
    real(real32), allocatable :: wf(:), hn(:), hr(:), hc(:)
    real(real64), allocatable :: hacc(:,:)
    integer(c_intptr_t) :: b_ptr, b_id, b_dist, b_w, b_bt6, b_n4, b_r4, kern
    integer(c_intptr_t) :: b_n4c, b_r4c
    real(kdkind) :: acc

    if (cfg%rank == 0) print *, 'Performing 4PCF (all configs, OpenCL bsearch)'

    nb     = cfg%nbins
    ncfg   = cfg%n_configs_4pcf
    n_hubs = cfg%num_data + cfg%num_rand
    call check_bt6_fits(nb)

    ngang = pick_ngang_4pcf(ncfg, 4)   ! 2 partials + 2 Kahan compensations
    ncol  = int(ncfg, int64) * ngang
    if (cfg%rank == 0) print '("  ngang=",i0,"  n_configs=",i0)', ngang, ncfg

    call pack_weights(wf)
    allocate(hn(ncol), hr(ncol), hc(ncol))
    hn = 0.0_real32
    hr = 0.0_real32
    hc = 0.0_real32
    allocate(hacc(ncol, 2))
    hacc = 0.0d0

    b_ptr  = cl_buf_in_i64(csr_ptr,  int(size(csr_ptr),int64))
    b_id   = cl_buf_in_i32(csr_id,   csr_total_edges)
    b_dist = cl_buf_in_i8 (csr_dist, csr_total_edges)
    b_w    = cl_buf_in_f32(wf,       n_hubs)
    b_bt6  = cl_buf_in_i32(bintable6, int(nb,int64)**6)   ! column-major == kernel flatten
    b_n4   = cl_buf_zeroed_f32(hn, ncol)
    b_r4   = cl_buf_zeroed_f32(hr, ncol)
    b_n4c  = cl_buf_zeroed_f32(hc, ncol)   ! Kahan compensation of b_n4
    b_r4c  = cl_buf_zeroed_f32(hc, ncol)   ! Kahan compensation of b_r4

    kern = cl_kernel_get('k_4pcf_all')
    call cl_arg_mem(kern, 0, b_ptr)
    call cl_arg_mem(kern, 1, b_id)
    call cl_arg_mem(kern, 2, b_dist)
    call cl_arg_mem(kern, 3, b_w)
    call cl_arg_mem(kern, 4, b_bt6)
    call cl_arg_mem(kern, 5, b_n4)
    call cl_arg_mem(kern, 6, b_r4)
    call cl_arg_i32(kern, 7,  nb)
    call cl_arg_i32(kern, 8,  ncfg)
    call cl_arg_i32(kern, 9, cfg%num_data)
    call cl_arg_i32(kern, 10, istart)
    call cl_arg_i32(kern, 11, iend)
    call cl_arg_i32(kern, 12, int(ngang, int32))
    call cl_arg_mem(kern, 17, b_n4c)
    call cl_arg_mem(kern, 18, b_r4c)
    if (cl_run_complete(kern, 13, 14, 15, 16, ngang)) then
      call cl_read_f32(b_n4, hn, ncol)
      call cl_read_f32(b_r4, hr, ncol)
      hacc(:, 1) = real(hn, real64)
      hacc(:, 2) = real(hr, real64)
      call cl_read_f32(b_n4c, hc, ncol)
      hacc(:, 1) = hacc(:, 1) - real(hc, real64)
      call cl_read_f32(b_r4c, hc, ncol)
      hacc(:, 2) = hacc(:, 2) - real(hc, real64)
    else
      if (cfg%rank == 0) print *, '  single GPU launch hit the watchdog; retrying tiled (slower)'
      call cl_write_f32(b_n4, hn, ncol)   ! hn/hr/hc are still all zeros here
      call cl_write_f32(b_r4, hr, ncol)
      call cl_write_f32(b_n4c, hc, ncol)
      call cl_write_f32(b_r4c, hc, ncol)
      call cl_run_bucketed(kern, 13, 14, 15, 16, ngang, n_hubs, &
                           [b_n4, b_n4c, b_r4, b_r4c], hacc)
    end if

    do c = 1, ncfg
      acc = 0.0d0
      do g = 0, int(ngang) - 1
        acc = acc + hacc(int(g,int64)*ncfg + c, 1)
      end do
      N4(c, 1) = N4(c, 1) + acc
      acc = 0.0d0
      do g = 0, int(ngang) - 1
        acc = acc + hacc(int(g,int64)*ncfg + c, 2)
      end do
      R4(c, 1) = R4(c, 1) + acc
    end do

    call cl_release(b_ptr); call cl_release(b_id); call cl_release(b_dist)
    call cl_release(b_w); call cl_release(b_bt6)
    call cl_release(b_n4); call cl_release(b_r4)
    call cl_release(b_n4c); call cl_release(b_r4c)
    call cl_release_kernel(kern)
    deallocate(wf, hn, hr, hc, hacc)

    call write_4pcf_results_noparity()
  end subroutine query_graph_4pcf_cl

  ! -------------------------------------------------------------------------
  ! 4PCF with parity decomposition.
  ! -------------------------------------------------------------------------
  subroutine query_graph_4pcf_parity_cl(istart, iend)
    integer, intent(in) :: istart, iend
    integer :: nb, ncfg, ndir, c, g
    integer :: p1, p2, p3
    integer(int64) :: ngang, n_hubs, ncol, idx
    real(real32), allocatable :: wf(:), hne(:), hno(:), hre(:), hro(:), hc(:)
    real(real64), allocatable :: hacc(:,:)
    integer(int8), allocatable :: signv(:)
    real(kdkind) :: vol, acc
    integer(c_intptr_t) :: b_ptr, b_id, b_dist, b_phi, b_w, b_bt6, b_sgn
    integer(c_intptr_t) :: b_ne, b_no, b_re, b_ro, kern
    integer(c_intptr_t) :: b_nec, b_noc, b_rec, b_roc

    if (cfg%rank == 0) print *, 'Performing 4PCF parity (all configs, OpenCL bsearch)'

    nb     = cfg%nbins
    ncfg   = cfg%n_configs_4pcf
    ndir   = cfg%n_dir_pixels
    n_hubs = cfg%num_data + cfg%num_rand
    call check_bt6_fits(nb)

    ngang = pick_ngang_4pcf(ncfg, 8)   ! 4 partials + 4 Kahan compensations
    ncol  = int(ncfg, int64) * ngang
    if (cfg%rank == 0) print '("  ngang=",i0,"  n_configs=",i0)', ngang, ncfg

    ! ---- Host-precompute the chirality sign table signv[(p1,p2,p3)] ---------
    ! Exactly the CPU formula in double, so the odd channel matches bit-for-bit.
    allocate(signv(int(ndir,int64)**3))
    do p1 = 1, ndir
      do p2 = 1, ndir
        do p3 = 1, ndir
          vol = dir_x(p1)*(dir_y(p2)*dir_z(p3) - dir_z(p2)*dir_y(p3)) &
              + dir_y(p1)*(dir_z(p2)*dir_x(p3) - dir_x(p2)*dir_z(p3)) &
              + dir_z(p1)*(dir_x(p2)*dir_y(p3) - dir_y(p2)*dir_x(p3))
          idx = int(p1-1,int64)*ndir*ndir + int(p2-1,int64)*ndir + (p3-1) + 1
          if (abs(vol) < VOL_DEGEN_TOL) then
            signv(idx) = 0_int8
          else if (vol > 0.0d0) then
            signv(idx) = 1_int8
          else
            signv(idx) = -1_int8
          end if
        end do
      end do
    end do

    call pack_weights(wf)
    allocate(hne(ncol), hno(ncol), hre(ncol), hro(ncol), hc(ncol))
    hne = 0.0_real32; hno = 0.0_real32; hre = 0.0_real32; hro = 0.0_real32
    hc = 0.0_real32
    allocate(hacc(ncol, 4))
    hacc = 0.0d0

    b_ptr  = cl_buf_in_i64(csr_ptr,  int(size(csr_ptr),int64))
    b_id   = cl_buf_in_i32(csr_id,   csr_total_edges)
    b_dist = cl_buf_in_i8 (csr_dist, csr_total_edges)
    b_phi  = cl_buf_in_i16(csr_phi,  csr_total_edges)
    b_w    = cl_buf_in_f32(wf,       n_hubs)
    b_bt6  = cl_buf_in_i32(bintable6, int(nb,int64)**6)
    b_sgn  = cl_buf_in_i8 (signv,    int(ndir,int64)**3)
    b_ne   = cl_buf_zeroed_f32(hne, ncol)
    b_no   = cl_buf_zeroed_f32(hno, ncol)
    b_re   = cl_buf_zeroed_f32(hre, ncol)
    b_ro   = cl_buf_zeroed_f32(hro, ncol)
    b_nec  = cl_buf_zeroed_f32(hc, ncol)   ! Kahan compensations of the four
    b_noc  = cl_buf_zeroed_f32(hc, ncol)   ! partials above, in the same order
    b_rec  = cl_buf_zeroed_f32(hc, ncol)
    b_roc  = cl_buf_zeroed_f32(hc, ncol)

    kern = cl_kernel_get('k_4pcf_parity')
    call cl_arg_mem(kern, 0, b_ptr)
    call cl_arg_mem(kern, 1, b_id)
    call cl_arg_mem(kern, 2, b_dist)
    call cl_arg_mem(kern, 3, b_phi)
    call cl_arg_mem(kern, 4, b_w)
    call cl_arg_mem(kern, 5, b_bt6)
    call cl_arg_mem(kern, 6, b_sgn)
    call cl_arg_mem(kern, 7,  b_ne)
    call cl_arg_mem(kern, 8,  b_no)
    call cl_arg_mem(kern, 9, b_re)
    call cl_arg_mem(kern, 10, b_ro)
    call cl_arg_i32(kern, 11, nb)
    call cl_arg_i32(kern, 12, ncfg)
    call cl_arg_i32(kern, 13, cfg%num_data)
    call cl_arg_i32(kern, 14, ndir)
    call cl_arg_i32(kern, 15, istart)
    call cl_arg_i32(kern, 16, iend)
    call cl_arg_i32(kern, 17, int(ngang, int32))
    call cl_arg_mem(kern, 22, b_nec)
    call cl_arg_mem(kern, 23, b_noc)
    call cl_arg_mem(kern, 24, b_rec)
    call cl_arg_mem(kern, 25, b_roc)
    if (cl_run_complete(kern, 18, 19, 20, 21, ngang)) then
      call cl_read_f32(b_ne, hne, ncol)
      call cl_read_f32(b_no, hno, ncol)
      call cl_read_f32(b_re, hre, ncol)
      call cl_read_f32(b_ro, hro, ncol)
      hacc(:, 1) = real(hne, real64)
      hacc(:, 2) = real(hno, real64)
      hacc(:, 3) = real(hre, real64)
      hacc(:, 4) = real(hro, real64)
      call cl_read_f32(b_nec, hc, ncol)
      hacc(:, 1) = hacc(:, 1) - real(hc, real64)
      call cl_read_f32(b_noc, hc, ncol)
      hacc(:, 2) = hacc(:, 2) - real(hc, real64)
      call cl_read_f32(b_rec, hc, ncol)
      hacc(:, 3) = hacc(:, 3) - real(hc, real64)
      call cl_read_f32(b_roc, hc, ncol)
      hacc(:, 4) = hacc(:, 4) - real(hc, real64)
    else
      if (cfg%rank == 0) print *, '  single GPU launch hit the watchdog; retrying tiled (slower)'
      call cl_write_f32(b_ne, hne, ncol)   ! host arrays are still all zeros here
      call cl_write_f32(b_no, hno, ncol)
      call cl_write_f32(b_re, hre, ncol)
      call cl_write_f32(b_ro, hro, ncol)
      call cl_write_f32(b_nec, hc, ncol)
      call cl_write_f32(b_noc, hc, ncol)
      call cl_write_f32(b_rec, hc, ncol)
      call cl_write_f32(b_roc, hc, ncol)
      call cl_run_bucketed(kern, 18, 19, 20, 21, ngang, n_hubs, &
                           [b_ne, b_nec, b_no, b_noc, b_re, b_rec, b_ro, b_roc], hacc)
    end if

    do c = 1, ncfg
      acc = 0.0d0
      do g = 0, int(ngang) - 1
        acc = acc + hacc(int(g,int64)*ncfg + c, 1)
      end do
      N4(c, 1) = N4(c, 1) + acc
      acc = 0.0d0
      do g = 0, int(ngang) - 1
        acc = acc + hacc(int(g,int64)*ncfg + c, 2)
      end do
      N4(c, 2) = N4(c, 2) + acc
      acc = 0.0d0
      do g = 0, int(ngang) - 1
        acc = acc + hacc(int(g,int64)*ncfg + c, 3)
      end do
      R4(c, 1) = R4(c, 1) + acc
      acc = 0.0d0
      do g = 0, int(ngang) - 1
        acc = acc + hacc(int(g,int64)*ncfg + c, 4)
      end do
      R4(c, 2) = R4(c, 2) + acc
    end do

    call cl_release(b_ptr); call cl_release(b_id); call cl_release(b_dist)
    call cl_release(b_phi); call cl_release(b_w)
    call cl_release(b_bt6); call cl_release(b_sgn)
    call cl_release(b_ne); call cl_release(b_no); call cl_release(b_re); call cl_release(b_ro)
    call cl_release(b_nec); call cl_release(b_noc); call cl_release(b_rec); call cl_release(b_roc)
    call cl_release_kernel(kern)
    deallocate(wf, hne, hno, hre, hro, hc, hacc, signv)

    call write_4pcf_results()
  end subroutine query_graph_4pcf_parity_cl

end module query_4pcf_cl_module
