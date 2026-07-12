! ===========================================================================
! GRAMSCI OpenCL driver — portable GPU build (Apple Silicon / Intel / AMD).
!
! Same structure and command-line interface as bin/gramsci and bin/gramsci_gpu:
! build the KD-tree graph on the CPU, flatten it to CSR, then offload the N-PCF
! counting to the GPU via OpenCL.  The 3PCF (isotropic), equilateral 3PCF, 4PCF
! and 4PCF-parity queries run on the GPU; 2PCF and any RSD (nmu>1) queries stay
! on the CPU, exactly as in the OpenACC build.
! ===========================================================================
program gramsci_cl
  use kdtree2_module
  use kdtree2_precision_module
  use node_module
  use config_module
  use io_module
  use graph_module
  use graph_utils_module
  use query_2pcf_module
  use query_3pcf_module
  use query_4pcf_module
  use csr_cl_module
  use cl_env_module
  use query_3pcf_cl_module
  use query_4pcf_cl_module
  implicit none
  INCLUDE "omp_lib.h"

  type(kdtree2), pointer :: kd_tree
  integer :: i, j, nn2, thread, threads
  real(kdkind) :: start, rand_val, avg_neighbors
  integer(8) :: wt0, wt1, wt_rate
  real(kdkind), allocatable :: sample_vec(:)

  thread = 0
  threads = 1

  ! ---- Configuration ----
  call default_params()
  cfg%backend = 'opencl'
  call parseOptions()

  ! ---- Initialise the OpenCL device (banner + kernel build) ----
  call cl_init()

  cfg%config_bins = cfg%nbins
  call create_binlookup()
  if (cfg%four_pcf_parity .or. cfg%four_pcf) call create_4pcf_binlookup()
  print *, 'Number of binned configurations:', cfg%config_bins
  print *, 'number of mu bins:', cfg%nmu

  if (cfg%rank == cfg%master) print *, 'allocated bins'

  ! ---- Count data ----
  call count_files_2()

  ! ---- Allocate data arrays ----
  allocate(points(cfg%d, cfg%num_data + cfg%num_rand))
  allocate(weights(cfg%num_data + cfg%num_rand))
  allocate(buffer(cfg%num_data + cfg%num_rand))
  allocate(radial_bins(cfg%nbins + 1))

  ! ---- Define radial bins ----
  print *, '================='
  print *, '= defining bins ='
  print *, '================='
  print *, cfg%rmin
  do i = 1, cfg%nbins + 1
    j = i - 1
    if (cfg%logbins) then
      radial_bins(i) = 10**(log10(cfg%rmin) + j * (log10(cfg%rmax) - log10(cfg%rmin)) / float(cfg%nbins))
    else
      radial_bins(i) = cfg%rmin + j * (cfg%rmax - cfg%rmin) / float(cfg%nbins)
    end if
    print *, radial_bins(i)
  end do

  ! ---- Read data ----
  call read_files_2()

  ! ---- Build KD-tree ----
  if (cfg%rank == 0) print *, 'building kd_tree '
  kd_tree => kdtree2_create(points, sort=.true., rearrange=.true.)

  if (cfg%rank == 0) print *, 'allocating arrays '
  call allocate_result_arrays()

  ! ---- Memory estimation ----
  allocate(sample_vec(100))
  do i = 1, 100
    call random_number(rand_val)
    j = floor(rand_val * cfg%num_data) + 1
    nn2 = kdtree2_r_count_around_point(tp=kd_tree, idxin=j, correltime=-1, r2=cfg%rmax*cfg%rmax)
    sample_vec(i) = real(nn2, kdkind)
  end do
  avg_neighbors = sum(sample_vec) / size(sample_vec)
  deallocate(sample_vec)

  ! ---- Compute bin parameters ----
  if (cfg%logbins) then
    cfg%delta_r = (log10(cfg%rmax) - log10(cfg%rmin)) / cfg%nbins
    cfg%inv_delta_r = 1.0d0 / cfg%delta_r
    cfg%log_rmin = log10(cfg%rmin)
  else
    cfg%delta_r = (cfg%rmax - cfg%rmin) / cfg%nbins
    cfg%inv_delta_r = 1.0d0 / cfg%delta_r
  end if
  cfg%mu_scale = 0.5d0 * cfg%nmu

  ! ---- Build graph ----
  if (cfg%rank == 0) print *, 'building relationships between nodes'
  allocate(output(cfg%num_data + cfg%num_rand))

  !$OMP PARALLEL
  !$ threads = OMP_GET_NUM_THREADS()
  !$ thread = OMP_GET_THREAD_NUM()
  !$OMP END PARALLEL

  if (cfg%rank == cfg%master) print *, 'Code running with ', threads, ' OMP threads'

  call set_kd_tree(kd_tree)

  call system_clock(wt0, wt_rate)
  call create_graph(1, cfg%num_data + cfg%num_rand)
  call system_clock(wt1)
  print '("Creating the graph took ",f12.3," seconds.")', &
    real(wt1 - wt0, kdkind) / real(wt_rate, kdkind)

  call kdtree2_destroy(kd_tree)
  deallocate(points)

  if (cfg%rank == 0) print *, 'finished building node relationships '

  ! ---- Phase 1: CPU queries that need the jagged output(:) ----
  call system_clock(wt0)
  N2 = 0.0d0
  N3 = 0.0d0

  ! Each query mode starts from freshly zeroed shared accumulators so that
  ! combined modes (-2pcf -3pcf ...) cannot cross-contaminate.

  ! 2PCF: CPU (O(N*m), fast; not worth GPU overhead)
  if (cfg%two_pcf) then
    N2 = 0.0d0
    N3 = 0.0d0
    call query_graph_2pcf(1, cfg%num_data + cfg%num_rand)
  end if

  ! Equilateral 3PCF with RSD: CPU only; isotropic runs on GPU below
  if (cfg%three_pcf_eq .and. cfg%RSD) then
    N2 = 0.0d0
    N3 = 0.0d0
    call query_graph_equilateral_triangle(1, cfg%num_data + cfg%num_rand)
  end if

  ! RSD (anisotropic) 3PCF: GPU version does not support nmu>1, run on CPU now
  if (cfg%three_pcf .and. cfg%RSD) then
    if (cfg%rank == 0) print *, 'RSD 3PCF: running on CPU'
    N2 = 0.0d0
    N3 = 0.0d0
    call query_graph_3pcf_all(1, cfg%num_data + cfg%num_rand)
  end if

  ! Internal 2PCF for disconnected 4PCF subtraction (and, in analytic -box
  ! mode, the 3PCF estimator subtraction): CPU, uses output(:)
  if (cfg%four_pcf .or. cfg%four_pcf_parity .or. &
      (cfg%analytic .and. (cfg%three_pcf .or. cfg%three_pcf_eq))) then
    call compute_2pcf_for_4pcf(1, cfg%num_data + cfg%num_rand)
  end if

  ! Analytic 4PCF: internal 3PCF pass (CPU, needs the jagged graph, so it
  ! must precede build_csr) storing zeta3_internal for the 4PCF estimator
  ! subtraction.  N2/N3 are cleared around the pass so the later OpenCL 3PCF
  ! queries start from clean accumulators.
  if (cfg%analytic .and. (cfg%four_pcf .or. cfg%four_pcf_parity)) then
    N2 = 0.0d0
    N3 = 0.0d0
    cfg%internal_3pcf = .true.
    call query_graph_3pcf_all(1, cfg%num_data + cfg%num_rand)
    cfg%internal_3pcf = .false.
    N2 = 0.0d0
    N3 = 0.0d0
  end if

  ! ---- Flatten graph to CSR (frees the jagged output(:)) ----
  call build_csr()
  if (cfg%rank == 0) print *, 'Jagged graph freed; using CSR from here'

  ! ---- Phase 2: GPU (OpenCL) queries using CSR ----
  if (cfg%three_pcf .and. .not. cfg%RSD) then
    N2 = 0.0d0
    N3 = 0.0d0
    call query_graph_3pcf_cl(1, cfg%num_data + cfg%num_rand)
  end if

  if (cfg%three_pcf_eq .and. .not. cfg%RSD) then
    N2 = 0.0d0
    N3 = 0.0d0
    call query_graph_equilateral_cl(1, cfg%num_data + cfg%num_rand)
  end if

  if (cfg%four_pcf) then
    allocate(N4(cfg%n_configs_4pcf, 1))
    allocate(R4(cfg%n_configs_4pcf, 1))
    N4 = 0.0d0 ; R4 = 0.0d0
    call query_graph_4pcf_cl(1, cfg%num_data + cfg%num_rand)
    deallocate(N4) ; deallocate(R4)
  end if

  if (cfg%four_pcf_parity) then
    call init_direction_lookup()
    allocate(N4(cfg%n_configs_4pcf, 2))
    allocate(R4(cfg%n_configs_4pcf, 2))
    N4 = 0.0d0 ; R4 = 0.0d0
    call query_graph_4pcf_parity_cl(1, cfg%num_data + cfg%num_rand)
    deallocate(N4) ; deallocate(R4)
    call cleanup_direction_lookup()
  end if

  call system_clock(wt1)
  print '("Querying graph took ",f12.3," seconds.")', &
    real(wt1 - wt0, kdkind) / real(wt_rate, kdkind)

  call deallocate_csr()
  call cl_shutdown()

  if (cfg%rank == 0) print *, 'finished querying the graph'

  call deallocate_arrays()
  deallocate(buffer)

  print *, "Exit... stage left!"

end program gramsci_cl
