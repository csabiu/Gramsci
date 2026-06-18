! ===========================================================================
! cl_module — minimal Fortran 2003 ISO_C_BINDING interface to the OpenCL 1.2 C
! API.  Hand-written (no third-party dependency); only the ~25 entry points and
! enums that the GRAMSCI OpenCL backend actually uses are declared.
!
! Opaque OpenCL handles (cl_platform_id, cl_device_id, cl_context,
! cl_command_queue, cl_program, cl_kernel, cl_mem, cl_event) are all
! pointer-sized and are represented here as integer(c_intptr_t).  On an LP64
! target (macOS arm64/x86_64, Linux x86_64) a c_intptr_t passed/returned by
! value occupies the same register as the C pointer, so this matches the ABI.
!
! Scalar cl_int/cl_uint are c_int32_t; cl_ulong / bitfields are c_int64_t;
! size_t is c_size_t.  void* arguments are passed as type(c_ptr) by value
! (obtained from c_loc of a host TARGET array/scalar).
!
! Tested against Apple's OpenCL 1.2 framework on Apple Silicon (M1) and the
! Khronos ICD loader on Linux.
! ===========================================================================
module cl_module
  use iso_c_binding
  implicit none
  public

  ! ---- Scalar type kinds (mnemonic aliases) ----
  integer, parameter :: cl_int    = c_int32_t
  integer, parameter :: cl_uint   = c_int32_t
  integer, parameter :: cl_long   = c_int64_t
  integer, parameter :: cl_ulong  = c_int64_t
  integer, parameter :: cl_bool   = c_int32_t
  integer, parameter :: cl_bitfield = c_int64_t

  ! ---- Status / boolean ----
  integer(c_int32_t), parameter :: CL_SUCCESS = 0
  integer(c_int32_t), parameter :: CL_TRUE  = 1
  integer(c_int32_t), parameter :: CL_FALSE = 0

  ! ---- cl_device_type (cl_bitfield) ----
  integer(c_int64_t), parameter :: CL_DEVICE_TYPE_DEFAULT     = 1_c_int64_t
  integer(c_int64_t), parameter :: CL_DEVICE_TYPE_CPU         = 2_c_int64_t
  integer(c_int64_t), parameter :: CL_DEVICE_TYPE_GPU         = 4_c_int64_t
  integer(c_int64_t), parameter :: CL_DEVICE_TYPE_ACCELERATOR = 8_c_int64_t
  integer(c_int64_t), parameter :: CL_DEVICE_TYPE_ALL         = int(z'FFFFFFFF', c_int64_t)

  ! ---- cl_mem_flags (cl_bitfield) ----
  integer(c_int64_t), parameter :: CL_MEM_READ_WRITE     = 1_c_int64_t       ! 1<<0
  integer(c_int64_t), parameter :: CL_MEM_WRITE_ONLY     = 2_c_int64_t       ! 1<<1
  integer(c_int64_t), parameter :: CL_MEM_READ_ONLY      = 4_c_int64_t       ! 1<<2
  integer(c_int64_t), parameter :: CL_MEM_USE_HOST_PTR   = 8_c_int64_t       ! 1<<3
  integer(c_int64_t), parameter :: CL_MEM_ALLOC_HOST_PTR = 16_c_int64_t      ! 1<<4
  integer(c_int64_t), parameter :: CL_MEM_COPY_HOST_PTR  = 32_c_int64_t      ! 1<<5

  ! ---- Platform info param names ----
  integer(c_int32_t), parameter :: CL_PLATFORM_NAME    = int(z'0902', c_int32_t)
  integer(c_int32_t), parameter :: CL_PLATFORM_VERSION = int(z'0901', c_int32_t)

  ! ---- Device info param names ----
  integer(c_int32_t), parameter :: CL_DEVICE_TYPE                = int(z'1000', c_int32_t)
  integer(c_int32_t), parameter :: CL_DEVICE_MAX_COMPUTE_UNITS   = int(z'1002', c_int32_t)
  integer(c_int32_t), parameter :: CL_DEVICE_MAX_WORK_GROUP_SIZE = int(z'1004', c_int32_t)
  integer(c_int32_t), parameter :: CL_DEVICE_MAX_MEM_ALLOC_SIZE  = int(z'1010', c_int32_t)
  integer(c_int32_t), parameter :: CL_DEVICE_GLOBAL_MEM_SIZE     = int(z'101F', c_int32_t)
  integer(c_int32_t), parameter :: CL_DEVICE_LOCAL_MEM_SIZE      = int(z'1023', c_int32_t)
  integer(c_int32_t), parameter :: CL_DEVICE_NAME                = int(z'102B', c_int32_t)
  integer(c_int32_t), parameter :: CL_DEVICE_DOUBLE_FP_CONFIG    = int(z'1032', c_int32_t)
  integer(c_int32_t), parameter :: CL_DEVICE_EXTENSIONS          = int(z'1030', c_int32_t)

  ! ---- Program build info param names ----
  integer(c_int32_t), parameter :: CL_PROGRAM_BUILD_STATUS = int(z'1181', c_int32_t)
  integer(c_int32_t), parameter :: CL_PROGRAM_BUILD_LOG    = int(z'1183', c_int32_t)

  interface

    ! cl_int clGetPlatformIDs(cl_uint, cl_platform_id*, cl_uint*)
    integer(c_int32_t) function clGetPlatformIDs(num_entries, platforms, num_platforms) &
        bind(c, name='clGetPlatformIDs')
      import :: c_int32_t, c_intptr_t
      integer(c_int32_t), value :: num_entries
      integer(c_intptr_t)       :: platforms(*)
      integer(c_int32_t)        :: num_platforms
    end function clGetPlatformIDs

    ! cl_int clGetPlatformInfo(cl_platform_id, cl_platform_info, size_t, void*, size_t*)
    integer(c_int32_t) function clGetPlatformInfo(platform, param_name, param_value_size, &
        param_value, param_value_size_ret) bind(c, name='clGetPlatformInfo')
      import :: c_int32_t, c_intptr_t, c_size_t, c_ptr
      integer(c_intptr_t), value :: platform
      integer(c_int32_t),  value :: param_name
      integer(c_size_t),   value :: param_value_size
      type(c_ptr),         value :: param_value
      integer(c_size_t)          :: param_value_size_ret
    end function clGetPlatformInfo

    ! cl_int clGetDeviceIDs(cl_platform_id, cl_device_type, cl_uint, cl_device_id*, cl_uint*)
    integer(c_int32_t) function clGetDeviceIDs(platform, device_type, num_entries, &
        devices, num_devices) bind(c, name='clGetDeviceIDs')
      import :: c_int32_t, c_int64_t, c_intptr_t
      integer(c_intptr_t), value :: platform
      integer(c_int64_t),  value :: device_type
      integer(c_int32_t),  value :: num_entries
      integer(c_intptr_t)        :: devices(*)
      integer(c_int32_t)         :: num_devices
    end function clGetDeviceIDs

    ! cl_int clGetDeviceInfo(cl_device_id, cl_device_info, size_t, void*, size_t*)
    integer(c_int32_t) function clGetDeviceInfo(device, param_name, param_value_size, &
        param_value, param_value_size_ret) bind(c, name='clGetDeviceInfo')
      import :: c_int32_t, c_intptr_t, c_size_t, c_ptr
      integer(c_intptr_t), value :: device
      integer(c_int32_t),  value :: param_name
      integer(c_size_t),   value :: param_value_size
      type(c_ptr),         value :: param_value
      integer(c_size_t)          :: param_value_size_ret
    end function clGetDeviceInfo

    ! cl_context clCreateContext(const cl_context_properties*, cl_uint, const cl_device_id*,
    !                            pfn_notify, void*, cl_int*)
    integer(c_intptr_t) function clCreateContext(properties, num_devices, devices, &
        pfn_notify, user_data, errcode_ret) bind(c, name='clCreateContext')
      import :: c_int32_t, c_intptr_t, c_ptr, c_funptr
      type(c_ptr),         value :: properties
      integer(c_int32_t),  value :: num_devices
      integer(c_intptr_t)        :: devices(*)
      type(c_funptr),      value :: pfn_notify
      type(c_ptr),         value :: user_data
      integer(c_int32_t)         :: errcode_ret
    end function clCreateContext

    ! cl_command_queue clCreateCommandQueue(cl_context, cl_device_id,
    !                                       cl_command_queue_properties, cl_int*)
    integer(c_intptr_t) function clCreateCommandQueue(context, device, properties, &
        errcode_ret) bind(c, name='clCreateCommandQueue')
      import :: c_int32_t, c_int64_t, c_intptr_t
      integer(c_intptr_t), value :: context
      integer(c_intptr_t), value :: device
      integer(c_int64_t),  value :: properties
      integer(c_int32_t)         :: errcode_ret
    end function clCreateCommandQueue

    ! cl_mem clCreateBuffer(cl_context, cl_mem_flags, size_t, void*, cl_int*)
    integer(c_intptr_t) function clCreateBuffer(context, flags, size, host_ptr, &
        errcode_ret) bind(c, name='clCreateBuffer')
      import :: c_int32_t, c_int64_t, c_intptr_t, c_size_t, c_ptr
      integer(c_intptr_t), value :: context
      integer(c_int64_t),  value :: flags
      integer(c_size_t),   value :: size
      type(c_ptr),         value :: host_ptr
      integer(c_int32_t)         :: errcode_ret
    end function clCreateBuffer

    ! cl_int clEnqueueWriteBuffer(cl_command_queue, cl_mem, cl_bool, size_t, size_t,
    !                             const void*, cl_uint, const cl_event*, cl_event*)
    integer(c_int32_t) function clEnqueueWriteBuffer(queue, buffer, blocking, offset, &
        size, ptr, num_events, event_wait_list, event) bind(c, name='clEnqueueWriteBuffer')
      import :: c_int32_t, c_intptr_t, c_size_t, c_ptr
      integer(c_intptr_t), value :: queue
      integer(c_intptr_t), value :: buffer
      integer(c_int32_t),  value :: blocking
      integer(c_size_t),   value :: offset
      integer(c_size_t),   value :: size
      type(c_ptr),         value :: ptr
      integer(c_int32_t),  value :: num_events
      type(c_ptr),         value :: event_wait_list
      type(c_ptr),         value :: event
    end function clEnqueueWriteBuffer

    ! cl_int clEnqueueReadBuffer(...)
    integer(c_int32_t) function clEnqueueReadBuffer(queue, buffer, blocking, offset, &
        size, ptr, num_events, event_wait_list, event) bind(c, name='clEnqueueReadBuffer')
      import :: c_int32_t, c_intptr_t, c_size_t, c_ptr
      integer(c_intptr_t), value :: queue
      integer(c_intptr_t), value :: buffer
      integer(c_int32_t),  value :: blocking
      integer(c_size_t),   value :: offset
      integer(c_size_t),   value :: size
      type(c_ptr),         value :: ptr
      integer(c_int32_t),  value :: num_events
      type(c_ptr),         value :: event_wait_list
      type(c_ptr),         value :: event
    end function clEnqueueReadBuffer

    ! cl_program clCreateProgramWithSource(cl_context, cl_uint, const char**,
    !                                      const size_t*, cl_int*)
    integer(c_intptr_t) function clCreateProgramWithSource(context, count, strings, &
        lengths, errcode_ret) bind(c, name='clCreateProgramWithSource')
      import :: c_int32_t, c_intptr_t, c_size_t, c_ptr
      integer(c_intptr_t), value :: context
      integer(c_int32_t),  value :: count
      type(c_ptr)                :: strings(*)
      integer(c_size_t)          :: lengths(*)
      integer(c_int32_t)         :: errcode_ret
    end function clCreateProgramWithSource

    ! cl_int clBuildProgram(cl_program, cl_uint, const cl_device_id*, const char*,
    !                       pfn_notify, void*)
    integer(c_int32_t) function clBuildProgram(program, num_devices, device_list, &
        options, pfn_notify, user_data) bind(c, name='clBuildProgram')
      import :: c_int32_t, c_intptr_t, c_char, c_funptr, c_ptr
      integer(c_intptr_t), value :: program
      integer(c_int32_t),  value :: num_devices
      integer(c_intptr_t)        :: device_list(*)
      character(kind=c_char)     :: options(*)
      type(c_funptr),      value :: pfn_notify
      type(c_ptr),         value :: user_data
    end function clBuildProgram

    ! cl_int clGetProgramBuildInfo(cl_program, cl_device_id, cl_program_build_info,
    !                              size_t, void*, size_t*)
    integer(c_int32_t) function clGetProgramBuildInfo(program, device, param_name, &
        param_value_size, param_value, param_value_size_ret) &
        bind(c, name='clGetProgramBuildInfo')
      import :: c_int32_t, c_intptr_t, c_size_t, c_ptr
      integer(c_intptr_t), value :: program
      integer(c_intptr_t), value :: device
      integer(c_int32_t),  value :: param_name
      integer(c_size_t),   value :: param_value_size
      type(c_ptr),         value :: param_value
      integer(c_size_t)          :: param_value_size_ret
    end function clGetProgramBuildInfo

    ! cl_kernel clCreateKernel(cl_program, const char*, cl_int*)
    integer(c_intptr_t) function clCreateKernel(program, kernel_name, errcode_ret) &
        bind(c, name='clCreateKernel')
      import :: c_int32_t, c_intptr_t, c_char
      integer(c_intptr_t), value :: program
      character(kind=c_char)     :: kernel_name(*)
      integer(c_int32_t)         :: errcode_ret
    end function clCreateKernel

    ! cl_int clSetKernelArg(cl_kernel, cl_uint, size_t, const void*)
    integer(c_int32_t) function clSetKernelArg(kernel, arg_index, arg_size, arg_value) &
        bind(c, name='clSetKernelArg')
      import :: c_int32_t, c_intptr_t, c_size_t, c_ptr
      integer(c_intptr_t), value :: kernel
      integer(c_int32_t),  value :: arg_index
      integer(c_size_t),   value :: arg_size
      type(c_ptr),         value :: arg_value
    end function clSetKernelArg

    ! cl_int clEnqueueNDRangeKernel(cl_command_queue, cl_kernel, cl_uint,
    !   const size_t*, const size_t*, const size_t*, cl_uint, const cl_event*, cl_event*)
    integer(c_int32_t) function clEnqueueNDRangeKernel(queue, kernel, work_dim, &
        global_work_offset, global_work_size, local_work_size, num_events, &
        event_wait_list, event) bind(c, name='clEnqueueNDRangeKernel')
      import :: c_int32_t, c_intptr_t, c_size_t, c_ptr
      integer(c_intptr_t), value :: queue
      integer(c_intptr_t), value :: kernel
      integer(c_int32_t),  value :: work_dim
      type(c_ptr),         value :: global_work_offset
      integer(c_size_t)          :: global_work_size(*)
      type(c_ptr),         value :: local_work_size
      integer(c_int32_t),  value :: num_events
      type(c_ptr),         value :: event_wait_list
      type(c_ptr),         value :: event
    end function clEnqueueNDRangeKernel

    integer(c_int32_t) function clFinish(queue) bind(c, name='clFinish')
      import :: c_int32_t, c_intptr_t
      integer(c_intptr_t), value :: queue
    end function clFinish

    integer(c_int32_t) function clReleaseMemObject(memobj) bind(c, name='clReleaseMemObject')
      import :: c_int32_t, c_intptr_t
      integer(c_intptr_t), value :: memobj
    end function clReleaseMemObject

    integer(c_int32_t) function clReleaseKernel(kernel) bind(c, name='clReleaseKernel')
      import :: c_int32_t, c_intptr_t
      integer(c_intptr_t), value :: kernel
    end function clReleaseKernel

    integer(c_int32_t) function clReleaseProgram(program) bind(c, name='clReleaseProgram')
      import :: c_int32_t, c_intptr_t
      integer(c_intptr_t), value :: program
    end function clReleaseProgram

    integer(c_int32_t) function clReleaseCommandQueue(queue) bind(c, name='clReleaseCommandQueue')
      import :: c_int32_t, c_intptr_t
      integer(c_intptr_t), value :: queue
    end function clReleaseCommandQueue

    integer(c_int32_t) function clReleaseContext(context) bind(c, name='clReleaseContext')
      import :: c_int32_t, c_intptr_t
      integer(c_intptr_t), value :: context
    end function clReleaseContext

  end interface

end module cl_module
