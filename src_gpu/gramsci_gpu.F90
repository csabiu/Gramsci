program Ngramsci_gpu
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
  use csr_module
  use query_3pcf_gpu_module
  use query_4pcf_gpu_module
  implicit none
  INCLUDE "omp_lib.h"

  type(kdtree2), pointer :: kd_tree
  integer :: i, j, thread, threads
  integer :: istart, iend
  real(kdkind) :: start, finish
  integer(8) :: wt0, wt1, wt_rate
  ! Wall-clock stage stamps ([t=...s] lines): the graph build and query
  ! have their own timers, but the untimed stages (catalogue parsing,
  ! KD-tree build/destroy, teardown) are where runtime-allocator overhead
  ! hides -- observed as an ~18 s mprotect storm under the nvfortran RTL.
  integer(8) :: t_stamp0, t_stamp_rate

  thread = 0
  threads = 1

  call system_clock(t_stamp0, t_stamp_rate)

  ! ---- Configuration ----
  call default_params()
  cfg%backend = 'openacc'
  call parseOptions()
  call stamp('options parsed')

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
  call stamp('catalogue count done; reading')
  call read_files_2()
  call stamp('catalogue read')

  ! ---- Build KD-tree ----
  if (cfg%rank == 0) print *, 'building kd_tree '
  kd_tree => kdtree2_create(points, sort=.true., rearrange=.true.)
  call stamp('kd-tree built')

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

  ! Graph build: time with system_clock (wall time).
  call system_clock(wt0, wt_rate)
  call cpu_time(start)
  call create_graph(1, 999)
  call cpu_time(finish)
  print '("Creating graph will take ~ ",f10.3," minutes.")', &
    (finish - start) * (cfg%num_data + cfg%num_rand) / (60.0 * 1000.0 * threads)
  call create_graph(1000, cfg%num_data + cfg%num_rand)
  call system_clock(wt1)
  print '("Creating the graph took ",f12.3," seconds.")', &
    real(wt1 - wt0, kdkind) / real(wt_rate, kdkind)

  call kdtree2_destroy(kd_tree)
  call stamp('kd-tree destroyed')
  if (cfg%exact_parity .and. cfg%four_pcf_parity) &
    call init_exact_parity_positions()
  deallocate(points)

  if (cfg%rank == 0) print *, 'finished building node relationships '

  ! ---- Phase 1: CPU queries that need the jagged output(:) ----
  ! These must all complete before build_csr() + output deallocation.
  call system_clock(wt0)
  N2 = 0.0d0
  N3 = 0.0d0

  ! 2PCF: CPU (O(N*m), fast; not worth GPU overhead)
  ! Each query mode starts from freshly zeroed shared accumulators so that
  ! combined modes (-2pcf -3pcf ...) cannot cross-contaminate.  N2jk/N3jk
  ! are shared across the pair/triplet modes, so they are re-zeroed too.
  if (cfg%two_pcf) then
    N2 = 0.0d0
    N3 = 0.0d0
    N2jk = 0.0d0
    N3jk = 0.0d0
    call query_graph_2pcf(1, cfg%num_data + cfg%num_rand)
  end if

  ! Equilateral 3PCF with RSD: CPU only (find_normal); isotropic runs on GPU
  ! below.  With -njk the isotropic case also runs here: the GPU equilateral
  ! kernel has no jackknife accumulation.
  if (cfg%three_pcf_eq .and. (cfg%RSD .or. cfg%njk > 0)) then
    if (cfg%njk > 0 .and. .not. cfg%RSD .and. cfg%rank == 0) &
      print *, 'equilateral 3PCF with -njk: running on CPU'
    N2 = 0.0d0
    N3 = 0.0d0
    N2jk = 0.0d0
    N3jk = 0.0d0
    call query_graph_equilateral_triangle(1, cfg%num_data + cfg%num_rand)
  end if

  ! RSD (anisotropic) 3PCF: GPU version does not support nmu>1, run on CPU now
  if (cfg%three_pcf .and. cfg%RSD) then
    if (cfg%rank == 0) print *, 'RSD 3PCF: running on CPU'
    N2 = 0.0d0
    N3 = 0.0d0
    N2jk = 0.0d0
    N3jk = 0.0d0
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
  ! subtraction.  N2/N3 are cleared around the pass so the later GPU 3PCF
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

  ! ---- Flatten graph to CSR ----
  ! build_csr frees each jagged row as it copies (and deallocates output),
  ! keeping peak host RAM at ~one graph copy.  All subsequent kernels use CSR.
  ! -merge runs no kernel (the raw sums come from the shard files), so the
  ! CSR is skipped and the jagged graph -- still needed above for the
  ! internal 2PCF -- is released here instead.
  if (cfg%merge_n > 0) then
    do i = 1, cfg%num_data + cfg%num_rand
      call output(i)%destroy()
    end do
    deallocate(output)
    if (cfg%rank == 0) print *, 'Jagged graph freed; -merge: no CSR, no GPU query'
  else
    call build_csr()
    if (cfg%rank == 0) print *, 'Jagged graph freed; using CSR from here'
  end if

  ! ---- Phase 2: GPU queries using CSR ----

  ! Isotropic 3PCF: GPU bsearch kernel (has slot-strided jackknife partials)
  if (cfg%three_pcf .and. .not. cfg%RSD) then
    N2 = 0.0d0
    N3 = 0.0d0
    N2jk = 0.0d0
    N3jk = 0.0d0
    call query_graph_3pcf_gpu(1, cfg%num_data + cfg%num_rand)
  end if

  ! Isotropic equilateral 3PCF: GPU kernel (RSD and -njk cases ran on CPU above)
  if (cfg%three_pcf_eq .and. .not. cfg%RSD .and. cfg%njk <= 0) then
    N2 = 0.0d0
    N3 = 0.0d0
    call query_graph_equilateral_gpu(1, cfg%num_data + cfg%num_rand)
  end if

  ! Both 4PCF kernels accumulate the jackknife touching sums on the device
  ! (direct atomics into N4jk/R4jk), so -njk no longer routes them to the
  ! CPU.  allocate_4pcf_jk is called unconditionally: the kernels' COPY
  ! clauses name N4jk/R4jk either way (dummy 1-region slice when off).
  if (cfg%four_pcf) then
    allocate(N4(cfg%n_configs_4pcf, 1))
    allocate(R4(cfg%n_configs_4pcf, 1))
    N4 = 0.0d0 ; R4 = 0.0d0
    call allocate_4pcf_jk(1)
    if (cfg%merge_n > 0) then
      call merge_4pcf_shards(.false.)
      call finish_4pcf_output(.false., 1, cfg%num_data + cfg%num_rand)
    else
      call query_hub_range(istart, iend)
      call query_graph_4pcf_gpu(istart, iend)
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
    if (cfg%merge_n > 0) then
      call merge_4pcf_shards(.true.)
      call finish_4pcf_output(.true., 1, cfg%num_data + cfg%num_rand)
    else
      call query_hub_range(istart, iend)
      call query_graph_4pcf_parity_gpu(istart, iend)
    end if
    deallocate(N4) ; deallocate(R4)
    call free_4pcf_jk()
    call cleanup_direction_lookup()
  end if

  ! Report wall-clock query time (system_clock captures GPU execution correctly;
  ! cpu_time/threads would show ~0 since the CPU idles while the GPU runs).
  call system_clock(wt1)
  print '("Querying graph took ",f12.3," seconds.")', &
    real(wt1 - wt0, kdkind) / real(wt_rate, kdkind)

  call deallocate_csr()

  if (cfg%rank == 0) print *, 'finished querying the graph'

  call deallocate_arrays()
  call stamp('teardown done')

  print *, "Exit... stage left!"

contains

  ! Hub range queried by this process: the whole catalogue, or under
  ! -shard k/n the k-th cost-balanced slice of the CSR rows.
  subroutine query_hub_range(istart, iend)
    integer, intent(out) :: istart, iend
    if (cfg%shard_n > 0) then
      call csr_shard_range(cfg%shard_k, cfg%shard_n, istart, iend)
    else
      istart = 1
      iend = cfg%num_data + cfg%num_rand
    end if
  end subroutine query_hub_range

  ! Elapsed wall time since program start, printed at stage boundaries.
  subroutine stamp(label)
    character(len=*), intent(in) :: label
    integer(8) :: t_now
    if (cfg%rank /= 0) return
    call system_clock(t_now)
    print '("[t=",f8.2," s] ",a)', &
      real(t_now - t_stamp0, kdkind) / real(t_stamp_rate, kdkind), label
  end subroutine stamp

end program Ngramsci_gpu
