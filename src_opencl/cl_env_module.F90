! ===========================================================================
! cl_env_module — OpenCL runtime environment for the GRAMSCI OpenCL backend.
!
! Owns the singleton platform/device/context/command-queue/program and exposes
! small, typed helpers so the query modules never touch raw ISO_C_BINDING:
!
!   call cl_init()                          ! pick a GPU, build the kernels
!   buf = cl_buf_in_i32(arr, n)             ! upload an array  -> cl_mem handle
!   buf = cl_buf_rw_f32(nbytes)             ! device-only RW scratch/output
!   k   = cl_kernel_get('name')             ! fetch a built kernel
!   call cl_arg_mem(k,0,buf) ; call cl_arg_i32(k,1,n) ; ...
!   call cl_run_1d(k, global_work)          ! enqueue + clFinish
!   call cl_read_f32(buf, host, n)          ! copy results back
!   call cl_release(buf)                    ! free a device buffer
!
! Device sizing facts queried at init are public module variables
! (cl_max_alloc, cl_global_mem, cl_has_fp64, cl_max_wg, cl_compute_units).
! ===========================================================================
module cl_env_module
  use iso_c_binding
  use iso_fortran_env, only: int8, int32, int64, real32, real64, error_unit
  use cl_module
  use cl_kernels_module, only: kernel_src
  implicit none
  private

  ! ---- Public handles & device facts ----
  integer(c_intptr_t), public :: cl_context_h = 0
  integer(c_intptr_t), public :: cl_queue_h   = 0
  integer(c_intptr_t), public :: cl_device_h  = 0
  integer(c_intptr_t), public :: cl_program_h = 0

  integer(c_int64_t),  public :: cl_max_alloc     = 0   ! bytes, largest single buffer
  integer(c_int64_t),  public :: cl_global_mem    = 0   ! bytes
  integer(c_int64_t),  public :: cl_max_wg        = 0   ! max work-group size
  integer(c_int32_t),  public :: cl_compute_units = 0
  logical,             public :: cl_has_fp64      = .false.
  character(len=256),  public :: cl_dev_name   = ''

  ! Number of interleaved hub buckets used by the tiled fallback (prime, so it
  ! is coprime to ngang).  Public so the kernel args and the launcher agree.
  integer, parameter, public :: CL_NBUCKETS = 1021

  public :: cl_init, cl_shutdown, cl_check
  public :: cl_buf_in_i8, cl_buf_in_i32, cl_buf_in_i64, cl_buf_in_f32
  public :: cl_buf_zeroed_f32, cl_buf_zeroed_i8
  public :: cl_read_f32, cl_read_i8, cl_write_f32, cl_release
  public :: cl_kernel_get, cl_release_kernel
  public :: cl_arg_mem, cl_arg_i32, cl_arg_i64, cl_arg_f32
  public :: cl_run_1d, cl_run_complete, cl_run_bucketed

contains

  ! -------------------------------------------------------------------------
  ! Initialise: choose the first GPU device (fall back to any device), create
  ! the context + in-order command queue, and build the embedded kernels.
  ! -------------------------------------------------------------------------
  subroutine cl_init()
    integer(c_int32_t) :: err, nplat, ndev
    integer(c_intptr_t) :: plats(16), devs(16)
    integer(c_intptr_t) :: plat
    integer :: ip
    integer(c_int64_t) :: fp64cfg
    logical :: found

    err = clGetPlatformIDs(16, plats, nplat)
    call cl_check(err, 'clGetPlatformIDs')
    if (nplat < 1) then
      write(error_unit,*) 'OpenCL: no platforms found'
      stop 1
    end if

    ! Prefer a platform that exposes a GPU; else take the first device anywhere.
    found = .false.
    plat = plats(1)
    do ip = 1, nplat
      err = clGetDeviceIDs(plats(ip), CL_DEVICE_TYPE_GPU, 16, devs, ndev)
      if (err == CL_SUCCESS .and. ndev >= 1) then
        plat = plats(ip)
        found = .true.
        exit
      end if
    end do
    if (.not. found) then
      do ip = 1, nplat
        err = clGetDeviceIDs(plats(ip), CL_DEVICE_TYPE_ALL, 16, devs, ndev)
        if (err == CL_SUCCESS .and. ndev >= 1) then
          plat = plats(ip)
          found = .true.
          exit
        end if
      end do
    end if
    if (.not. found) then
      write(error_unit,*) 'OpenCL: no devices found on any platform'
      stop 1
    end if
    cl_device_h = devs(1)

    ! ---- Query device facts ----
    cl_dev_name   = dev_info_str(cl_device_h, CL_DEVICE_NAME)
    cl_max_alloc     = dev_info_i64(cl_device_h, CL_DEVICE_MAX_MEM_ALLOC_SIZE)
    cl_global_mem    = dev_info_i64(cl_device_h, CL_DEVICE_GLOBAL_MEM_SIZE)
    cl_max_wg        = dev_info_size(cl_device_h, CL_DEVICE_MAX_WORK_GROUP_SIZE)
    cl_compute_units = dev_info_u32(cl_device_h, CL_DEVICE_MAX_COMPUTE_UNITS)
    fp64cfg          = dev_info_i64(cl_device_h, CL_DEVICE_DOUBLE_FP_CONFIG)
    cl_has_fp64      = (fp64cfg /= 0_c_int64_t)

    print '(a)',        ' OpenCL device: '//trim(cl_dev_name)
    print '(a,i0,a,f6.2,a)', '   compute units: ', cl_compute_units, &
          '   max buffer: ', real(cl_max_alloc)/1.0e9, ' GB'
    print '(a,f6.2,a,l1)', '   global mem: ', real(cl_global_mem)/1.0e9, &
          ' GB   fp64 in kernels: ', cl_has_fp64

    ! ---- Context + command queue ----
    cl_context_h = clCreateContext(c_null_ptr, 1, devs, c_null_funptr, c_null_ptr, err)
    call cl_check(err, 'clCreateContext')
    cl_queue_h = clCreateCommandQueue(cl_context_h, cl_device_h, 0_c_int64_t, err)
    call cl_check(err, 'clCreateCommandQueue')

    call build_program()
  end subroutine cl_init

  ! -------------------------------------------------------------------------
  ! Build the embedded OpenCL C program; on failure dump the build log.
  ! -------------------------------------------------------------------------
  subroutine build_program()
    integer(c_int32_t) :: err
    character(kind=c_char), allocatable, target :: src(:)
    character(kind=c_char), target :: opts(1)
    type(c_ptr) :: strings(1)
    integer(c_size_t) :: lengths(1)
    integer(c_intptr_t) :: devs(1)
    character(len=:), allocatable :: ksrc
    integer :: i, n

    ksrc = kernel_src()
    n = len(ksrc)
    allocate(src(n))
    do i = 1, n
      src(i) = ksrc(i:i)
    end do
    strings(1) = c_loc(src(1))
    lengths(1) = int(n, c_size_t)

    cl_program_h = clCreateProgramWithSource(cl_context_h, 1, strings, lengths, err)
    call cl_check(err, 'clCreateProgramWithSource')

    opts(1) = c_null_char
    devs(1) = cl_device_h
    err = clBuildProgram(cl_program_h, 1, devs, opts, c_null_funptr, c_null_ptr)
    if (err /= CL_SUCCESS) then
      write(error_unit,'(a,i0,a)') 'OpenCL: clBuildProgram failed (', err, '). Build log:'
      call dump_build_log()
      stop 1
    end if
    deallocate(src)
  end subroutine build_program

  subroutine dump_build_log()
    integer(c_int32_t) :: err
    integer(c_size_t) :: logsz, ret
    character(kind=c_char), allocatable, target :: logbuf(:)
    character(len=:), allocatable :: fstr
    integer :: i

    err = clGetProgramBuildInfo(cl_program_h, cl_device_h, CL_PROGRAM_BUILD_LOG, &
                                0_c_size_t, c_null_ptr, logsz)
    if (logsz <= 1) return
    allocate(logbuf(logsz))
    err = clGetProgramBuildInfo(cl_program_h, cl_device_h, CL_PROGRAM_BUILD_LOG, &
                                logsz, c_loc(logbuf(1)), ret)
    allocate(character(len=int(logsz)) :: fstr)
    do i = 1, int(logsz)
      fstr(i:i) = logbuf(i)
    end do
    write(error_unit,'(a)') trim(fstr)
    deallocate(logbuf, fstr)
  end subroutine dump_build_log

  ! -------------------------------------------------------------------------
  ! Error checking / decoding
  ! -------------------------------------------------------------------------
  subroutine cl_check(err, where)
    integer(c_int32_t), intent(in) :: err
    character(*), intent(in) :: where
    if (err == CL_SUCCESS) return
    write(error_unit,'(a,a,a,i0,a,a)') 'OpenCL ERROR in ', where, ' : code ', err, &
                                       '  ', trim(cl_err_str(err))
    stop 1
  end subroutine cl_check

  pure function cl_err_str(e) result(s)
    integer(c_int32_t), intent(in) :: e
    character(len=40) :: s
    select case (e)
    case (-1);  s = 'DEVICE_NOT_FOUND'
    case (-2);  s = 'DEVICE_NOT_AVAILABLE'
    case (-3);  s = 'COMPILER_NOT_AVAILABLE'
    case (-4);  s = 'MEM_OBJECT_ALLOCATION_FAILURE'
    case (-5);  s = 'OUT_OF_RESOURCES'
    case (-6);  s = 'OUT_OF_HOST_MEMORY'
    case (-11); s = 'BUILD_PROGRAM_FAILURE'
    case (-30); s = 'INVALID_VALUE'
    case (-34); s = 'INVALID_CONTEXT'
    case (-36); s = 'INVALID_COMMAND_QUEUE'
    case (-38); s = 'INVALID_MEM_OBJECT'
    case (-44); s = 'INVALID_PROGRAM'
    case (-45); s = 'INVALID_PROGRAM_EXECUTABLE'
    case (-46); s = 'INVALID_KERNEL_NAME'
    case (-48); s = 'INVALID_KERNEL'
    case (-49); s = 'INVALID_ARG_INDEX'
    case (-50); s = 'INVALID_ARG_VALUE'
    case (-51); s = 'INVALID_ARG_SIZE'
    case (-52); s = 'INVALID_KERNEL_ARGS'
    case (-54); s = 'INVALID_WORK_GROUP_SIZE'
    case (-55); s = 'INVALID_WORK_ITEM_SIZE'
    case (-61); s = 'INVALID_BUFFER_SIZE'
    case default; s = 'UNKNOWN'
    end select
  end function cl_err_str

  ! -------------------------------------------------------------------------
  ! Buffer creation — upload host arrays (READ_ONLY | COPY_HOST_PTR)
  ! -------------------------------------------------------------------------
  function cl_buf_in_i8(arr, n) result(buf)
    integer(int8), target, intent(in) :: arr(*)
    integer(int64), intent(in) :: n
    integer(c_intptr_t) :: buf
    integer(c_int32_t) :: err
    buf = clCreateBuffer(cl_context_h, ior(CL_MEM_READ_ONLY, CL_MEM_COPY_HOST_PTR), &
                         int(n, c_size_t), c_loc(arr(1)), err)
    call cl_check(err, 'clCreateBuffer(i8)')
  end function cl_buf_in_i8

  function cl_buf_in_i32(arr, n) result(buf)
    integer(int32), target, intent(in) :: arr(*)
    integer(int64), intent(in) :: n
    integer(c_intptr_t) :: buf
    integer(c_int32_t) :: err
    buf = clCreateBuffer(cl_context_h, ior(CL_MEM_READ_ONLY, CL_MEM_COPY_HOST_PTR), &
                         int(n*4_int64, c_size_t), c_loc(arr(1)), err)
    call cl_check(err, 'clCreateBuffer(i32)')
  end function cl_buf_in_i32

  function cl_buf_in_i64(arr, n) result(buf)
    integer(int64), target, intent(in) :: arr(*)
    integer(int64), intent(in) :: n
    integer(c_intptr_t) :: buf
    integer(c_int32_t) :: err
    buf = clCreateBuffer(cl_context_h, ior(CL_MEM_READ_ONLY, CL_MEM_COPY_HOST_PTR), &
                         int(n*8_int64, c_size_t), c_loc(arr(1)), err)
    call cl_check(err, 'clCreateBuffer(i64)')
  end function cl_buf_in_i64

  function cl_buf_in_f32(arr, n) result(buf)
    real(real32), target, intent(in) :: arr(*)
    integer(int64), intent(in) :: n
    integer(c_intptr_t) :: buf
    integer(c_int32_t) :: err
    buf = clCreateBuffer(cl_context_h, ior(CL_MEM_READ_ONLY, CL_MEM_COPY_HOST_PTR), &
                         int(n*4_int64, c_size_t), c_loc(arr(1)), err)
    call cl_check(err, 'clCreateBuffer(f32)')
  end function cl_buf_in_f32

  ! Device read/write float buffer initialised from a zeroed host array (so the
  ! accumulators start at exactly 0.0 without a separate fill kernel).
  function cl_buf_zeroed_f32(arr, n) result(buf)
    real(real32), target, intent(in) :: arr(*)
    integer(int64), intent(in) :: n
    integer(c_intptr_t) :: buf
    integer(c_int32_t) :: err
    buf = clCreateBuffer(cl_context_h, ior(CL_MEM_READ_WRITE, CL_MEM_COPY_HOST_PTR), &
                         int(n*4_int64, c_size_t), c_loc(arr(1)), err)
    call cl_check(err, 'clCreateBuffer(zeroed_f32)')
  end function cl_buf_zeroed_f32

  subroutine cl_read_f32(buf, arr, n)
    integer(c_intptr_t), intent(in) :: buf
    real(real32), target, intent(inout) :: arr(*)
    integer(int64), intent(in) :: n
    integer(c_int32_t) :: err
    err = clEnqueueReadBuffer(cl_queue_h, buf, CL_TRUE, 0_c_size_t, &
                              int(n*4_int64, c_size_t), c_loc(arr(1)), 0, &
                              c_null_ptr, c_null_ptr)
    call cl_check(err, 'clEnqueueReadBuffer(f32)')
  end subroutine cl_read_f32

  ! Overwrite a device float buffer from a host array (used to re-zero
  ! accumulators before the tiled fallback re-runs).
  subroutine cl_write_f32(buf, arr, n)
    integer(c_intptr_t), intent(in) :: buf
    real(real32), target, intent(in) :: arr(*)
    integer(int64), intent(in) :: n
    integer(c_int32_t) :: err
    err = clEnqueueWriteBuffer(cl_queue_h, buf, CL_TRUE, 0_c_size_t, &
                               int(n*4_int64, c_size_t), c_loc(arr(1)), 0, &
                               c_null_ptr, c_null_ptr)
    call cl_check(err, 'clEnqueueWriteBuffer(f32)')
  end subroutine cl_write_f32

  ! Device read/write int8 buffer initialised from a host array (used for the
  ! per-work-item completion flags).
  function cl_buf_zeroed_i8(arr, n) result(buf)
    integer(int8), target, intent(in) :: arr(*)
    integer(int64), intent(in) :: n
    integer(c_intptr_t) :: buf
    integer(c_int32_t) :: err
    buf = clCreateBuffer(cl_context_h, ior(CL_MEM_READ_WRITE, CL_MEM_COPY_HOST_PTR), &
                         int(n, c_size_t), c_loc(arr(1)), err)
    call cl_check(err, 'clCreateBuffer(zeroed_i8)')
  end function cl_buf_zeroed_i8

  subroutine cl_read_i8(buf, arr, n)
    integer(c_intptr_t), intent(in) :: buf
    integer(int8), target, intent(inout) :: arr(*)
    integer(int64), intent(in) :: n
    integer(c_int32_t) :: err
    err = clEnqueueReadBuffer(cl_queue_h, buf, CL_TRUE, 0_c_size_t, &
                              int(n, c_size_t), c_loc(arr(1)), 0, &
                              c_null_ptr, c_null_ptr)
    call cl_check(err, 'clEnqueueReadBuffer(i8)')
  end subroutine cl_read_i8

  subroutine cl_release(buf)
    integer(c_intptr_t), intent(inout) :: buf
    integer(c_int32_t) :: err
    if (buf == 0) return
    err = clReleaseMemObject(buf)
    call cl_check(err, 'clReleaseMemObject')
    buf = 0
  end subroutine cl_release

  ! -------------------------------------------------------------------------
  ! Kernels & arguments
  ! -------------------------------------------------------------------------
  function cl_kernel_get(name) result(k)
    character(*), intent(in) :: name
    integer(c_intptr_t) :: k
    integer(c_int32_t) :: err
    character(kind=c_char), allocatable, target :: cname(:)
    integer :: i, n
    n = len_trim(name)
    allocate(cname(n+1))
    do i = 1, n
      cname(i) = name(i:i)
    end do
    cname(n+1) = c_null_char
    k = clCreateKernel(cl_program_h, cname, err)
    call cl_check(err, 'clCreateKernel('//trim(name)//')')
    deallocate(cname)
  end function cl_kernel_get

  subroutine cl_release_kernel(kernel)
    integer(c_intptr_t), intent(inout) :: kernel
    integer(c_int32_t) :: err
    if (kernel == 0) return
    err = clReleaseKernel(kernel)
    kernel = 0
  end subroutine cl_release_kernel

  subroutine cl_arg_mem(kernel, idx, mem)
    integer(c_intptr_t), intent(in) :: kernel
    integer, intent(in) :: idx
    integer(c_intptr_t), intent(in) :: mem
    integer(c_intptr_t), target :: m
    integer(c_int32_t) :: err
    m = mem
    err = clSetKernelArg(kernel, int(idx, c_int32_t), &
                         int(c_sizeof(m), c_size_t), c_loc(m))
    call cl_check(err, 'clSetKernelArg(mem)')
  end subroutine cl_arg_mem

  subroutine cl_arg_i32(kernel, idx, ival)
    integer(c_intptr_t), intent(in) :: kernel
    integer, intent(in) :: idx
    integer(int32), intent(in) :: ival
    integer(c_int32_t), target :: v
    integer(c_int32_t) :: err
    v = ival
    err = clSetKernelArg(kernel, int(idx, c_int32_t), 4_c_size_t, c_loc(v))
    call cl_check(err, 'clSetKernelArg(i32)')
  end subroutine cl_arg_i32

  subroutine cl_arg_i64(kernel, idx, lval)
    integer(c_intptr_t), intent(in) :: kernel
    integer, intent(in) :: idx
    integer(int64), intent(in) :: lval
    integer(c_int64_t), target :: v
    integer(c_int32_t) :: err
    v = lval
    err = clSetKernelArg(kernel, int(idx, c_int32_t), 8_c_size_t, c_loc(v))
    call cl_check(err, 'clSetKernelArg(i64)')
  end subroutine cl_arg_i64

  subroutine cl_arg_f32(kernel, idx, fval)
    integer(c_intptr_t), intent(in) :: kernel
    integer, intent(in) :: idx
    real(real32), intent(in) :: fval
    real(real32), target :: v
    integer(c_int32_t) :: err
    v = fval
    err = clSetKernelArg(kernel, int(idx, c_int32_t), 4_c_size_t, c_loc(v))
    call cl_check(err, 'clSetKernelArg(f32)')
  end subroutine cl_arg_f32

  ! -------------------------------------------------------------------------
  ! Launch a 1-D NDRange (let the runtime choose the local work-group size)
  ! and block until it finishes.
  ! -------------------------------------------------------------------------
  subroutine cl_run_1d(kernel, global)
    integer(c_intptr_t), intent(in) :: kernel
    integer(int64), intent(in) :: global
    integer(c_size_t) :: gws(1)
    integer(c_int32_t) :: err
    gws(1) = int(global, c_size_t)
    err = clEnqueueNDRangeKernel(cl_queue_h, kernel, 1, c_null_ptr, gws, &
                                 c_null_ptr, 0, c_null_ptr, c_null_ptr)
    call cl_check(err, 'clEnqueueNDRangeKernel')
    err = clFinish(cl_queue_h)
    call cl_check(err, 'clFinish')
  end subroutine cl_run_1d

  ! Fast path: run the kernel as ONE launch over ALL hubs (bucket window
  ! [0,NBUCKETS) -> the kernel's full-window fast path skips the bucket modulo).
  ! Each work-item sets done_flag[g]=1 only after sweeping all its hubs, so if
  ! the GPU watchdog truncates the launch (Apple's OpenCL then returns partial
  ! results with NO error) some flags stay 0.  Returns .true. iff every
  ! work-item finished.  A single sustained launch is by far the fastest mode on
  ! Apple's GPU: many short, clFinish-synced launches keep the GPU clock from
  ! ramping up (DVFS), so fragmenting can be several times slower.  Hence we try
  ! one launch first and only fall back to tiling (cl_run_bucketed) if it was
  ! truncated.  Arg indices are 0-based.
  function cl_run_complete(kernel, ia_nbuckets, ia_blo, ia_bhi, ia_flag, global) &
      result(complete)
    integer(c_intptr_t), intent(in) :: kernel
    integer, intent(in) :: ia_nbuckets, ia_blo, ia_bhi, ia_flag
    integer(int64), intent(in) :: global
    logical :: complete
    integer(int8), allocatable, target :: flags(:)
    integer(c_intptr_t) :: bflag

    allocate(flags(global))
    flags = 0_int8
    bflag = cl_buf_zeroed_i8(flags, global)
    call cl_arg_mem(kernel, ia_flag, bflag)
    call cl_arg_i32(kernel, ia_nbuckets, CL_NBUCKETS)
    call cl_arg_i32(kernel, ia_blo, 0)
    call cl_arg_i32(kernel, ia_bhi, CL_NBUCKETS)
    call cl_run_1d(kernel, global)
    call cl_read_i8(bflag, flags, global)
    complete = all(flags == 1_int8)
    call cl_release(bflag)
    deallocate(flags)
  end function cl_run_complete

  ! Fallback: run the kernel over all hubs as a sequence of INTERLEAVED bucket
  ! windows, each launch short enough to clear the GPU watchdog.  Used only when
  ! cl_run_complete reported truncation (heavy 4PCF, etc.); the caller must
  ! re-zero the accumulators first (the truncated launch left partial counts).
  !
  ! The split is by *bucket* (hub i -> bucket i mod NBUCKETS), NOT by contiguous
  ! hub range: per-hub cost is very non-uniform (clustered data hubs cost ~50x
  ! the uniform random hubs), so a contiguous slice would isolate all the
  ! expensive hubs into one launch and under-occupy the GPU.  Interleaved
  ! buckets give every launch a balanced, catalog-wide spread.  NBUCKETS is
  ! prime, hence coprime to ngang (a power of two), so each work-item's hubs
  ! span many buckets and partial windows stay well occupied.
  !
  ! Window sizing is adaptive: seed the first window near ~80% occupancy from
  ! the hubs/work-item ratio, then time each launch and grow toward TARGET_T
  ! seconds (capped per step so a launch cannot overshoot the watchdog).
  ! Override with GRAMSCI_CL_TARGET_SEC.
  subroutine cl_run_bucketed(kernel, ia_nbuckets, ia_blo, ia_bhi, ia_flag, global, n_hubs)
    integer(c_intptr_t), intent(in) :: kernel
    integer, intent(in) :: ia_nbuckets, ia_blo, ia_bhi, ia_flag
    integer(int64), intent(in) :: global, n_hubs
    integer :: blo, bhi, gb, doneb, gb0
    integer(int64) :: t0, t1, rate
    real(real64) :: dt, target_t, hpw
    integer(int8), allocatable, target :: flags(:)
    integer(c_intptr_t) :: bflag
    character(32) :: env
    integer :: stat

    target_t = 0.35d0
    call get_environment_variable('GRAMSCI_CL_TARGET_SEC', env, status=stat)
    if (stat == 0) read(env, *) target_t
    if (target_t <= 0.0d0) target_t = 0.35d0

    ! flag buffer is required by the kernel signature; not inspected here since
    ! adaptively-sized windows stay under the watchdog.
    allocate(flags(global))
    flags = 0_int8
    bflag = cl_buf_zeroed_i8(flags, global)
    call cl_arg_mem(kernel, ia_flag, bflag)
    call cl_arg_i32(kernel, ia_nbuckets, CL_NBUCKETS)

    hpw = real(n_hubs, real64) / real(max(global, 1_int64), real64)
    if (hpw < 1.0d0) hpw = 1.0d0
    gb0 = nint(CL_NBUCKETS * (1.0d0 - 0.2d0 ** (1.0d0 / hpw)))
    gb0 = max(min(gb0, 96), 1)

    blo = 0
    gb  = gb0
    do while (blo < CL_NBUCKETS)
      bhi = min(blo + gb, CL_NBUCKETS)
      call cl_arg_i32(kernel, ia_blo, blo)
      call cl_arg_i32(kernel, ia_bhi, bhi)

      call system_clock(t0, rate)
      call cl_run_1d(kernel, global)
      call system_clock(t1)
      dt    = real(t1 - t0, real64) / real(rate, real64)
      doneb = bhi - blo

      if (dt > 1.0d-4) then
        gb = int(real(doneb, real64) * target_t / dt)  ! buckets to fill target
        gb = min(gb, doneb * 4)                         ! cap growth per step
        gb = max(gb, 1)
      else
        gb = gb * 4                                     ! too fast to time; ramp
      end if
      gb = min(gb, CL_NBUCKETS)

      blo = bhi
    end do

    call cl_release(bflag)
    deallocate(flags)
  end subroutine cl_run_bucketed

  ! -------------------------------------------------------------------------
  ! Teardown
  ! -------------------------------------------------------------------------
  subroutine cl_shutdown()
    integer(c_int32_t) :: err
    if (cl_program_h /= 0) err = clReleaseProgram(cl_program_h)
    if (cl_queue_h   /= 0) err = clReleaseCommandQueue(cl_queue_h)
    if (cl_context_h /= 0) err = clReleaseContext(cl_context_h)
    cl_program_h = 0; cl_queue_h = 0; cl_context_h = 0
  end subroutine cl_shutdown

  ! -------------------------------------------------------------------------
  ! Device-info query helpers
  ! -------------------------------------------------------------------------
  function dev_info_i64(dev, param) result(v)
    integer(c_intptr_t), intent(in) :: dev
    integer(c_int32_t), intent(in) :: param
    integer(c_int64_t), target :: v
    integer(c_int32_t) :: err
    integer(c_size_t) :: ret
    v = 0
    err = clGetDeviceInfo(dev, param, 8_c_size_t, c_loc(v), ret)
    call cl_check(err, 'clGetDeviceInfo(i64)')
  end function dev_info_i64

  function dev_info_size(dev, param) result(v)
    integer(c_intptr_t), intent(in) :: dev
    integer(c_int32_t), intent(in) :: param
    integer(c_size_t), target :: tmp
    integer(c_int64_t) :: v
    integer(c_int32_t) :: err
    integer(c_size_t) :: ret
    tmp = 0
    err = clGetDeviceInfo(dev, param, int(c_sizeof(tmp), c_size_t), c_loc(tmp), ret)
    call cl_check(err, 'clGetDeviceInfo(size_t)')
    v = int(tmp, c_int64_t)
  end function dev_info_size

  function dev_info_u32(dev, param) result(v)
    integer(c_intptr_t), intent(in) :: dev
    integer(c_int32_t), intent(in) :: param
    integer(c_int32_t), target :: v
    integer(c_int32_t) :: err
    integer(c_size_t) :: ret
    v = 0
    err = clGetDeviceInfo(dev, param, 4_c_size_t, c_loc(v), ret)
    call cl_check(err, 'clGetDeviceInfo(u32)')
  end function dev_info_u32

  function dev_info_str(dev, param) result(s)
    integer(c_intptr_t), intent(in) :: dev
    integer(c_int32_t), intent(in) :: param
    character(len=256) :: s
    character(kind=c_char), target :: buf(256)
    integer(c_int32_t) :: err
    integer(c_size_t) :: ret
    integer :: i
    s = ''
    buf = c_null_char
    err = clGetDeviceInfo(dev, param, 256_c_size_t, c_loc(buf(1)), ret)
    if (err /= CL_SUCCESS) return
    do i = 1, min(256, int(ret))
      if (buf(i) == c_null_char) exit
      s(i:i) = buf(i)
    end do
  end function dev_info_str

end module cl_env_module
