program Ngramsci
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
  implicit none
  INCLUDE "omp_lib.h"

  type(kdtree2), pointer :: kd_tree
  integer :: i, j, thread, threads
  real(kdkind) :: start, finish
  logical :: use_bsearch

  thread = 0
  threads = 1

  ! ---- Configuration ----
  call default_params()
  call parseOptions()

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
  call read_jk_regions()
  call allocate_result_arrays()

  ! ---- Memory estimation (first-principles; samples the whole catalogue) ----
  call print_graph_ram_estimate(kd_tree)

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

  ! cfg%half_graph = .true. is the default; only edges with neighbor_id > hub_id are stored.
  call cpu_time(start)
  call create_graph(1, 999)
  call cpu_time(finish)
  print '("Creating graph will take ~ ",f10.3," minutes.")', &
    (finish - start) * (cfg%num_data + cfg%num_rand) / (60.0 * 1000.0 * threads)
  print *, 'If this takes longer than time to drink a coffee, maybe you should give me more CPUs'
  call create_graph(1000, cfg%num_data + cfg%num_rand)
  call cpu_time(finish)
  print '("Creating the graph took ",f12.3," seconds.")', (finish - start) / threads

  call kdtree2_destroy(kd_tree)
  if (cfg%exact_parity .and. cfg%four_pcf_parity) &
    call init_exact_parity_positions()
  deallocate(points)

  if (cfg%rank == 0) print *, 'finished building node relationships '

  ! ---- Query the graph ----
  call cpu_time(start)
  N2 = 0.0d0
  N3 = 0.0d0

  ! GRAMSCI_BSEARCH=1 selects the legacy binary-search kernels (benchmarking
  ! the v1 inner loop against the merge-walk; not exposed on the CLI).
  block
    character(8) :: env_bs
    integer :: env_stat
    call get_environment_variable('GRAMSCI_BSEARCH', env_bs, status=env_stat)
    use_bsearch = (env_stat == 0)
    if (use_bsearch .and. cfg%rank == 0) &
      print *, 'NOTE: using legacy binary-search kernels (GRAMSCI_BSEARCH set)'
  end block

  ! Each query mode starts from freshly zeroed shared accumulators so that
  ! combined modes (-2pcf -3pcf ...) cannot cross-contaminate.  The
  ! jackknife touching sums share N2jk/N3jk across the pair/triplet modes,
  ! so they are re-zeroed with them.
  if (cfg%two_pcf) then
    N2 = 0.0d0
    N3 = 0.0d0
    N2jk = 0.0d0
    N3jk = 0.0d0
    call query_graph_2pcf(1, cfg%num_data + cfg%num_rand)
  end if

  ! Internal 2PCF: needed for the 4PCF disconnected term and, in analytic
  ! (-box, no -ran) mode, for the 3PCF/equilateral estimator subtraction.
  ! Must run before the 3PCF queries, whose write routines use xi_2pcf.
  if (cfg%four_pcf .or. cfg%four_pcf_parity .or. &
      (cfg%analytic .and. (cfg%three_pcf .or. cfg%three_pcf_eq))) then
    call compute_2pcf_for_4pcf(1, cfg%num_data + cfg%num_rand)
  end if

  if (cfg%three_pcf) then
    N2 = 0.0d0
    N3 = 0.0d0
    N2jk = 0.0d0
    N3jk = 0.0d0
    if (use_bsearch) then
      call query_graph_3pcf_all_bsearch(1, cfg%num_data + cfg%num_rand)
    else
      call query_graph_3pcf_all(1, cfg%num_data + cfg%num_rand)
    end if
  end if
  if (cfg%three_pcf_eq) then
    N2 = 0.0d0
    N3 = 0.0d0
    N2jk = 0.0d0
    N3jk = 0.0d0
    call query_graph_equilateral_triangle(1, cfg%num_data + cfg%num_rand)
  end if

  ! Analytic 4PCF: internal 3PCF pass storing the connected zeta3 per config
  ! (zeta3_internal) for the 4PCF estimator subtraction.  N2/N3 are cleared
  ! around the pass so it neither inherits nor leaks triangle counts.
  if (cfg%analytic .and. (cfg%four_pcf .or. cfg%four_pcf_parity)) then
    N2 = 0.0d0
    N3 = 0.0d0
    cfg%internal_3pcf = .true.
    call query_graph_3pcf_all(1, cfg%num_data + cfg%num_rand)
    cfg%internal_3pcf = .false.
    N2 = 0.0d0
    N3 = 0.0d0
  end if

  if (cfg%four_pcf) then
    allocate(N4(cfg%n_configs_4pcf, 1))
    allocate(R4(cfg%n_configs_4pcf, 1))
    N4 = 0.0d0 ; R4 = 0.0d0
    call allocate_4pcf_jk(1)
    if (use_bsearch) then
      call query_graph_4pcf_bsearch(1, cfg%num_data + cfg%num_rand)
    else
      call query_graph_4pcf(1, cfg%num_data + cfg%num_rand)
    end if
    deallocate(N4) ; deallocate(R4)
    call free_4pcf_jk()
  end if

  if (cfg%four_pcf_parity) then
    call init_direction_lookup()
    allocate(N4(cfg%n_configs_4pcf, 2))
    allocate(R4(cfg%n_configs_4pcf, 2))
    N4 = 0.0d0 ; R4 = 0.0d0
    call allocate_4pcf_jk(2)
    if (use_bsearch) then
      call query_graph_4pcf_parity_bsearch(1, cfg%num_data + cfg%num_rand)
    else
      call query_graph_4pcf_parity(1, cfg%num_data + cfg%num_rand)
    end if
    deallocate(N4) ; deallocate(R4)
    call free_4pcf_jk()
    call cleanup_direction_lookup()
  end if

  ! ---- Cleanup graph ----
  do i = 1, size(output)
    call output(i)%destroy()
  end do
  deallocate(output)

  call cpu_time(finish)
  print '("Querying graph took ",f12.3," seconds.")', (finish - start) / threads

  if (cfg%rank == 0) print *, 'finished querying the graph'

  call deallocate_arrays()


  print *, "Exit... stage left!"

end program Ngramsci
