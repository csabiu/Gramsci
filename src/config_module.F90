module config_module
  use kdtree2_precision_module
  use node_module
  use iso_fortran_env, only: int8
  implicit none

  ! Single source of truth for the release version (see CHANGELOG.md);
  ! printed by -version and stamped into every output file's header.
  character(len=*), parameter :: GRAMSCI_VERSION = '2.5.0'

  ! Configuration derived type bundling all parameters
  type :: gramsci_config
    integer :: d = 3
    integer :: nbins = 0
    integer :: nmu = 1
    integer :: config_bins = 0
    real(kdkind) :: rmin = 0.0d0
    real(kdkind) :: rmax = 0.0d0
    real(kdkind) :: delta_r, inv_delta_r, log_rmin
    real(kdkind) :: mu_scale
    logical :: logbins = .false.
    logical :: RSD = .false.
    logical :: rancat = .false.
    logical :: two_pcf = .false.
    logical :: three_pcf = .false.
    logical :: three_pcf_eq = .false.
    logical :: four_pcf = .false.
    logical :: four_pcf_parity = .false.
    logical :: exact_parity = .false.
    logical :: half_graph = .true.
    integer :: n_configs_4pcf = 0
    ! Direction pixel parameters for 4PCF parity (n_theta * n_phi <= 32767:
    ! the pixel index is stored as int16 per graph edge).  Override with
    ! -ntheta / -nphi; ignored entirely under -exactparity.
    !
    ! The default 8 x 32 = 256 is the coarsest grid at which the recovered
    ! parity-odd amplitude is unbiased for every tetrahedron shape: on the
    ! historical 4 x 16 grid the near-degenerate (flattened, sheared) shapes
    ! were attenuated by ~13% and ~9%.  The finer grid costs ~1% of query
    ! time and one extra byte per graph edge.
    integer :: n_theta_dir = 8
    integer :: n_phi_dir = 32
    integer :: n_dir_pixels = 256   ! = n_theta_dir * n_phi_dir
    ! Periodic-box mode (-box): minimum-image separations; if no random
    ! catalogue is given, RR/RRR/RRRR are computed analytically.
    logical :: periodic = .false.
    logical :: analytic = .false.
    ! When true, the 3PCF write routine stores zeta3_internal (per config)
    ! for the 4PCF connected-part subtraction instead of writing a file.
    logical :: internal_3pcf = .false.
    ! More than one query mode requested: each mode writes to
    ! <output_file>.<mode> instead of sharing output_file.
    logical :: multi_mode = .false.
    real(kdkind) :: boxsize(3) = 0.0d0
    ! Power sums of the normalized data weights (set in read_files_2),
    ! used for the analytic pair/triple/quadruple weight normalizations.
    real(kdkind) :: sw2 = 0.0d0, sw3 = 0.0d0, sw4 = 0.0d0
    ! Accumulate pair Legendre multipoles (L2, L4 of the pair mu) during
    ! graph construction, for the anisotropic disconnected-4PCF subtraction.
    ! In a periodic box the line of sight is plane-parallel (z axis); in
    ! survey mode it is the pair midpoint direction (requires -nmu > 1).
    logical :: disc_rsd = .false.
    character(len=2000) :: file1 = '', file2 = ''
    character(len=2000) :: output_file = 'result.txt'
    ! Which build is running ('cpu' | 'openacc' | 'opencl'); set by the
    ! driver before parseOptions, reported by -version and the file headers
    character(len=8) :: backend = 'cpu'
    ! Process identity: always rank 0 of 1 (kept so the write routines'
    ! rank guards read uniformly across the CPU/GPU/OpenCL drivers)
    integer :: rank = 0
    integer :: num_tasks = 1
    integer :: master = 0
    ! Data dimensions
    integer :: num_data = 0
    integer :: num_rand = 0
    ! Delete-one jackknife. Region labels are supplied per catalogue, 1..njk.
    ! A triplet contributes to realisation m unless one of its three points
    ! lies in region m, so a single pass yields all njk realisations:
    !   N_m = N_total - (sum over triplets touching region m)
    integer :: njk = 0
    character(len=2000) :: jkgal = '', jkran = ''
    ! No -jkgal/-jkran given with -njk: partition the sky internally into
    ! njk equal-count angular regions (see assign_jk_regions_angular).
    logical :: jk_internal = .false.
  end type gramsci_config

  ! Singleton config instance
  type(gramsci_config), save :: cfg

  ! Shared data arrays (accessible to all modules)
  type(node), dimension(:), allocatable :: output
  real(kdkind), dimension(:,:), allocatable :: points
  real(kdkind), allocatable :: weights(:), radial_bins(:)
  real(kdkind), allocatable :: N2(:,:,:), N3(:,:,:)
  ! Per-region sums of triplets TOUCHING each region (not the realisations
  ! themselves); the realisations are formed by subtraction at write time.
  real(kdkind), allocatable :: N2jk(:,:,:,:), N3jk(:,:,:,:)
  integer, allocatable :: region(:)
  integer, allocatable :: bintable(:,:,:,:)
  ! 4PCF parity arrays
  real(kdkind), allocatable :: N4(:,:), R4(:,:)
  ! Per-region sums of quadruplets TOUCHING each region, shape
  ! (n_configs_4pcf, nch, njk) with nch = 1 (-4pcf) or 2 (-4pcfp).
  ! Updated with !$OMP ATOMIC, not a reduction (see allocate_4pcf_jk).
  real(kdkind), allocatable :: N4jk(:,:,:), R4jk(:,:,:)
  integer, allocatable :: bintable6(:,:,:,:,:,:)
  ! Canonical bin 6-tuple for each 4PCF config index (shape [6, n_configs_4pcf])
  integer, allocatable :: canon_bins_4pcf(:,:)
  ! Number of distinct ordered 6-tuples in each 4PCF config's S4 orbit
  ! (filled by create_4pcf_binlookup; needed by the analytic RRRR)
  integer, allocatable :: orbit_mult_4pcf(:)
  ! Internal 2PCF for disconnected 4PCF subtraction
  real(kdkind), allocatable :: DD_2pcf(:), RR_2pcf(:), xi_2pcf(:)
  ! Per-region touching sums of the internal 2PCF pair counts (nbins, njk),
  ! used to form the per-realisation xi0 in the jackknifed zeta_disc.
  real(kdkind), allocatable :: DD_2pcf_jk(:,:), RR_2pcf_jk(:,:)
  ! Pair Legendre-multipole sums (accumulated in create_graph when
  ! cfg%disc_rsd) and the resulting xi_2 / xi_4 multipoles, used for the
  ! anisotropic disconnected-4PCF subtraction in redshift space.
  real(kdkind), allocatable :: sum_pair_l2(:), sum_pair_l4(:)
  real(kdkind), allocatable :: xi2_2pcf(:), xi4_2pcf(:)
  ! Internal connected 3PCF per config (analytic mode, 4PCF subtraction)
  real(kdkind), allocatable :: zeta3_internal(:)

contains

  subroutine default_params()
    cfg%d = 3
    cfg%logbins = .false.
    cfg%nbins = 0
    cfg%rmin = 0.0d0
    cfg%rmax = 0.0d0
    cfg%output_file = 'result.txt'
    cfg%RSD = .false.
    cfg%nmu = 1
    cfg%rancat = .false.
    cfg%three_pcf_eq = .false.
    cfg%three_pcf = .false.
    cfg%four_pcf = .false.
    cfg%four_pcf_parity = .false.
    cfg%two_pcf = .false.
    cfg%half_graph = .true.
    cfg%periodic = .false.
    cfg%analytic = .false.
    cfg%internal_3pcf = .false.
    cfg%boxsize = 0.0d0
    cfg%sw2 = 0.0d0
    cfg%sw3 = 0.0d0
    cfg%sw4 = 0.0d0
    cfg%disc_rsd = .false.
  end subroutine default_params

  subroutine parseOptions()
    use extension, only: getArgument, nArguments
    implicit none
    integer(kind=4) :: i, n
    integer :: ios
    character(len=32) :: opt
    character(len=2000) :: arg

    n = nArguments()

    ! -version must work standalone, so scan for it before any arg-count check
    do i = 1, n
      call getArgument(i, opt)
      if (opt == '-version') then
        print '(a)', 'gramsci ' // GRAMSCI_VERSION // ' (' // trim(cfg%backend) // ' build)'
        stop
      end if
    end do

    if (n < 6 .and. cfg%rank == cfg%master) then
      print *, ' '
      print *, 'Not enough input parameters. Please read the following help info'
      print *, ' '
      call print_help()
      stop 1
    end if

    i = 1
    do while (i <= n)
      call getArgument(i, opt)
      select case (opt)
      case ('-gal')
        call get_value(arg)
        cfg%file1 = trim(arg)
        i = i + 2
      case ('-ran')
        cfg%rancat = .true.
        call get_value(arg)
        cfg%file2 = trim(arg)
        i = i + 2
      case ('-jkgal')
        call getArgument(i+1, arg)
        cfg%jkgal = trim(arg)
        i = i + 2
      case ('-jkran')
        call getArgument(i+1, arg)
        cfg%jkran = trim(arg)
        i = i + 2
      case ('-njk')
        call getArgument(i+1, arg)
        read(arg, *) cfg%njk
        i = i + 2
      case ('-out')
        call get_value(arg)
        cfg%output_file = trim(arg)
        i = i + 2
      case ('-rmin')
        call get_value(arg)
        read(arg, *, iostat=ios) cfg%rmin
        call check_numeric()
        i = i + 2
      case ('-rmax')
        call get_value(arg)
        read(arg, *, iostat=ios) cfg%rmax
        call check_numeric()
        i = i + 2
      case ('-nbins')
        call get_value(arg)
        read(arg, *, iostat=ios) cfg%nbins
        call check_numeric()
        i = i + 2
      case ('-nmu')
        call get_value(arg)
        read(arg, *, iostat=ios) cfg%nmu
        call check_numeric()
        if (cfg%nmu >= 2) then
          cfg%RSD = .true.
          if (cfg%rank == cfg%master) print *, 'Anisotropic analysis requested'
        end if
        i = i + 2
      case ('-box')
        call get_value(arg)
        call parse_boxsize(trim(arg))
        cfg%periodic = .true.
        i = i + 2
      case ('-log')
        cfg%logbins = .true.
        if (cfg%rank == cfg%master) print *, 'Using logarithmically spaced bins'
        i = i + 1
      case ('-help')
        call print_help()
        stop
      case ('-version')   ! handled by the pre-scan above; never reached
        i = i + 1
      case ('-3pcf')
        cfg%three_pcf = .true.
        i = i + 1
      case ('-4pcf')
        cfg%four_pcf = .true.
        i = i + 1
      case ('-equi')
        cfg%three_pcf_eq = .true.
        i = i + 1
      case ('-2pcf')
        cfg%two_pcf = .true.
        i = i + 1
      case ('-4pcfp')
        cfg%four_pcf_parity = .true.
        i = i + 1
      case ('-exactparity')
        cfg%exact_parity = .true.
        i = i + 1
      case ('-ntheta')
        call get_value(arg)
        read(arg, *, iostat=ios) cfg%n_theta_dir
        call check_numeric()
        i = i + 2
      case ('-nphi')
        call get_value(arg)
        read(arg, *, iostat=ios) cfg%n_phi_dir
        call check_numeric()
        i = i + 2
      case default
        print '("ERROR: unknown option ",a)', trim(opt)
        stop 1
      end select
    end do

    ! --- Post-parse validation ---
    ! Query modes may be combined (the graph is built once and each query
    ! starts from freshly zeroed accumulators).  With more than one mode the
    ! outputs go to <output_file>.<mode> so they cannot clobber each other.
    block
      integer :: n_modes
      n_modes = count([cfg%two_pcf, cfg%three_pcf, cfg%three_pcf_eq, &
                       cfg%four_pcf, cfg%four_pcf_parity])
      if (n_modes == 0) then
        print *, 'ERROR: specify at least one query mode (-2pcf | -3pcf | -equi | -4pcf | -4pcfp)'
        stop 1
      end if
      cfg%multi_mode = n_modes > 1
      if (cfg%multi_mode .and. cfg%rank == cfg%master) then
        print *, 'Multiple query modes: writing one output file per mode', &
                 ' (' // trim(cfg%output_file) // '.<mode>)'
      end if
    end block
    if (cfg%nbins <= 0) then
      print *, 'ERROR: -nbins must be specified and > 0'
      stop 1
    end if
    if (cfg%rmax <= cfg%rmin) then
      print *, 'ERROR: -rmax must be greater than -rmin'
      stop 1
    end if
    if (cfg%logbins .and. cfg%rmin <= 0.0d0) then
      print *, 'ERROR: -log requires -rmin > 0'
      stop 1
    end if
    if (cfg%nbins > 127) then
      print *, 'ERROR: -nbins must be <= 127 (bin indices stored as int8)'
      stop 1
    end if
    if (cfg%nmu > 127) then
      print *, 'ERROR: -nmu must be <= 127 (mu indices stored as int8)'
      stop 1
    end if
    if (cfg%n_theta_dir < 1 .or. cfg%n_phi_dir < 1) then
      print *, 'ERROR: -ntheta and -nphi must be >= 1'
      stop 1
    end if
    ! Division form: the product would overflow default integer for large
    ! -ntheta/-nphi and could wrap to a small positive that passes the test.
    if (cfg%n_theta_dir > 32767 / cfg%n_phi_dir) then
      print *, 'ERROR: -ntheta * -nphi must be <= 32767 (phi pixel stored as int16)'
      stop 1
    end if
    ! Derived AFTER validation: sole sizing input for dir_x/dir_y/dir_z.
    cfg%n_dir_pixels = cfg%n_theta_dir * cfg%n_phi_dir

    ! --- Periodic-box mode validation ---
    if (cfg%periodic) then
      if (any(cfg%boxsize <= 0.0d0)) then
        print *, 'ERROR: -box requires positive box side length(s)'
        stop 1
      end if
      if (cfg%nmu > 1) then
        print *, 'ERROR: -box supports isotropic analysis only (-nmu 1);'
        print *, '       the midpoint line-of-sight mu is not defined in a periodic box'
        stop 1
      end if
      if (cfg%rmax >= 0.5d0 * minval(cfg%boxsize)) then
        print *, 'ERROR: -box requires rmax < L/2 (minimum-image uniqueness)'
        stop 1
      end if
      if ((cfg%three_pcf .or. cfg%three_pcf_eq .or. cfg%four_pcf .or. cfg%four_pcf_parity) &
          .and. cfg%rmax > 0.25d0 * minval(cfg%boxsize)) then
        print *, 'ERROR: 3PCF/4PCF with -box require rmax <= L/4 so that the'
        print *, '       minimum-image side lengths of every tuple are mutually consistent'
        stop 1
      end if
      cfg%analytic = .not. cfg%rancat
      if (cfg%rank == cfg%master) then
        print '("Periodic box: ",3(f12.4,1x))', cfg%boxsize
        if (cfg%analytic) then
          print *, 'No random catalogue given: RR counts will be computed analytically'
        else
          print *, 'Random catalogue given: using catalogue randoms with periodic distances'
        end if
      end if
    end if

    ! --- Jackknife validation ---
    if (cfg%njk < 0 .or. cfg%njk == 1) then
      print *, 'ERROR: -njk must be >= 2 (delete-one needs at least 2 regions)'
      stop 1
    end if
    if (cfg%njk > 0) then
      cfg%jk_internal = (len_trim(cfg%jkgal) == 0 .and. len_trim(cfg%jkran) == 0)
      if (cfg%analytic) then
        print *, 'ERROR: -njk is incompatible with analytic randoms (-box without'
        print *, '       -ran): the delete-one RR/RRR/RRRR of a region is undefined'
        stop 1
      end if
      if (.not. cfg%rancat) then
        print *, 'ERROR: -njk requires -ran (delete-one random counts)'
        stop 1
      end if
      if (cfg%periodic .and. cfg%jk_internal) then
        print *, 'ERROR: internal angular jackknife regions are undefined in a'
        print *, '       periodic box (no observer direction); supply -jkgal/-jkran'
        stop 1
      end if
      if (cfg%jk_internal .and. cfg%rank == cfg%master) then
        print '(" jackknife: partitioning the sky internally into ",i0, &
               &" equal-count angular regions")', cfg%njk
      end if
    end if

    ! Anisotropic disconnected-4PCF subtraction: needs a per-pair mu, which
    ! exists in a periodic box (plane-parallel z line of sight) or in survey
    ! mode when RSD (-nmu > 1) is on.  Costs a few percent at graph build.
    cfg%disc_rsd = (cfg%four_pcf .or. cfg%four_pcf_parity) .and. &
                   (cfg%periodic .or. cfg%RSD)
    if (cfg%disc_rsd .and. cfg%rank == cfg%master) then
      print *, '4PCF: accumulating pair multipoles (xi_2, xi_4) for the'
      print *, '      anisotropic disconnected-term subtraction'
    end if

  contains

    ! Fetch the value of option opt (argument i+1), erroring out if absent.
    subroutine get_value(val)
      character(len=*), intent(out) :: val
      if (i + 1 > n) then
        print '("ERROR: option ",a," requires a value")', trim(opt)
        stop 1
      end if
      call getArgument(i + 1, val)
    end subroutine get_value

    subroutine check_numeric()
      if (ios /= 0) then
        print '("ERROR: cannot parse value for option ",a,": ",a)', trim(opt), trim(arg)
        stop 1
      end if
    end subroutine check_numeric


  end subroutine parseOptions

  ! Parse the -box argument: either a single side length L (cubic box)
  ! or three comma-separated values Lx,Ly,Lz.
  subroutine parse_boxsize(arg)
    character(len=*), intent(in) :: arg
    integer :: c1, c2, ios
    c1 = index(arg, ',')
    if (c1 == 0) then
      read(arg, *, iostat=ios) cfg%boxsize(1)
      cfg%boxsize(2:3) = cfg%boxsize(1)
    else
      c2 = c1 + index(arg(c1+1:), ',')
      if (c2 == c1) then
        print *, 'ERROR: -box takes either L or Lx,Ly,Lz'
        stop 1
      end if
      read(arg(:c1-1), *, iostat=ios) cfg%boxsize(1)
      if (ios == 0) read(arg(c1+1:c2-1), *, iostat=ios) cfg%boxsize(2)
      if (ios == 0) read(arg(c2+1:), *, iostat=ios) cfg%boxsize(3)
    end if
    if (ios /= 0) then
      print *, 'ERROR: cannot parse -box argument: ', trim(arg)
      stop 1
    end if
  end subroutine parse_boxsize

  subroutine print_help()
    print *, 'PURPOSE: Code for calculating the N-PCF of a 3D point set.'
    print *, ' '
    print *, 'CALLING SEQUENCE:'
    print *, '      gramsci [-gal galaxy file] [-ran ranfile (optional)]'
    print *, '              [-rmin Rmin] [-rmax Rmax] [-nbins Nbins]'
    print *, '              [-nmu Nmu] [-out out_file]'
    print *, '              [-2pcf | -3pcf | -equi | -4pcf | -4pcfp]'
    print *, '              [-log] [-box L]'
    print *, ' '
    print *, '      eg: gramsci -gal test.gal -ran test.ran -rmin 10.0 -rmax 30.0 -nbins 10 -2pcf'
    print *, ' '
    print *, 'INPUTS:'
    print *, '       -gal   galaxy catalogue file (x y z [weight])'
    print *, '       -ran   random catalogue file (same format)'
    print *, '       -out   output filename'
    print *, '       -rmin  min radial separation'
    print *, '       -rmax  max radial separation'
    print *, ' '
    print *, 'OPTIONAL:'
    print *, '       -nbins  number of radial bins'
    print *, '       -nmu    number of mu bins (enables anisotropic analysis)'
    print *, '       -log    use logarithmic radial binning'
    print *, '       -version  print the version and exit'
    print *, '       -box    periodic box side length L (or Lx,Ly,Lz).'
    print *, '               Uses minimum-image separations; without -ran the'
    print *, '               RR/RRR/RRRR counts are computed analytically'
    print *, '       -exactparity  parity sign from the exact galaxy positions'
    print *, '               instead of pixelized spoke directions (-4pcfp only)'
    print *, '       -ntheta N     polar direction bins   (default 8,  -4pcfp only)'
    print *, '       -nphi   M     azimuthal direction bins (default 32, N*M <= 32767)'
    print *, ' '
    print *, 'JACKKNIFE ERRORS (delete-one, all query modes; requires -ran):'
    print *, '       -njk N  number of jackknife regions.  Without -jkgal/-jkran'
    print *, '               the sky (as seen from the origin) is partitioned'
    print *, '               internally into N equal-count angular regions'
    print *, '               (sin(dec) bands x phi slices; boundaries from the'
    print *, '               randoms).  Purely angular by construction, so radial'
    print *, '               and angular systematics are never mixed.  The labels'
    print *, '               used are written to <out>.jkgal / <out>.jkran.'
    print *, '       -jkgal  file with one region label (1..N) per galaxy row'
    print *, '       -jkran  file with one region label (1..N) per random row'
    print *, '       Each mode writes <out>[.<mode>].jk (realisations) and'
    print *, '       <out>[.<mode>].jkerr (jackknife mean and sigma per bin).'
    print *, ' '
    print *, 'QUERY MODES (combinable; the graph is built once per run):'
    print *, '       -2pcf   2-point correlation function'
    print *, '       -3pcf   3-point correlation function (all configs)'
    print *, '       -equi   3-point correlation function (equilateral only)'
    print *, '       -4pcf   4-point correlation function (all configs)'
    print *, '       -4pcfp  4-point correlation function (all configs + parity)'
    print *, '       With one mode the result goes to -out exactly; with'
    print *, '       several, each mode writes to <out>.<mode> (e.g. res.3pcf)'
    print *, ' '
  end subroutine print_help

  subroutine create_binlookup()
    integer :: i, j, k

    ! The analytic 4PCF estimator needs the internal per-config 3PCF for its
    ! connected-part subtraction, so build the 3PCF lookup in that mode too.
    if (cfg%three_pcf .or. &
        (cfg%analytic .and. (cfg%four_pcf .or. cfg%four_pcf_parity))) then
      allocate(bintable(cfg%nbins, cfg%nbins, cfg%nbins, 1))
      cfg%config_bins = 0
      do i = 1, cfg%nbins
        do j = i, cfg%nbins
          do k = j, cfg%nbins
            cfg%config_bins = cfg%config_bins + 1
            bintable(i, j, k, 1) = cfg%config_bins
            bintable(i, k, j, 1) = cfg%config_bins
            bintable(j, i, k, 1) = cfg%config_bins
            bintable(j, k, i, 1) = cfg%config_bins
            bintable(k, i, j, 1) = cfg%config_bins
            bintable(k, j, i, 1) = cfg%config_bins
          end do
        end do
      end do
    end if
  end subroutine create_binlookup

  subroutine allocate_result_arrays()
    implicit none
    allocate(N2(cfg%config_bins, cfg%nmu, 3))
    allocate(N3(cfg%config_bins, cfg%nmu, 3))
    N2 = 0.0d0
    N3 = 0.0d0
    ! Always allocate: the OpenMP reduction clause in the query loop names
    ! these arrays unconditionally, and reducing over an unallocated
    ! allocatable segfaults. When jackknife is off they are a single dummy
    ! slice that is never written (the guards are on cfg%njk > 0).
    allocate(N2jk(cfg%config_bins, cfg%nmu, 3, max(cfg%njk, 1)))
    allocate(N3jk(cfg%config_bins, cfg%nmu, 3, max(cfg%njk, 1)))
    N2jk = 0.0d0
    N3jk = 0.0d0
    if (cfg%njk > 0 .and. cfg%rank == 0) &
      print *, 'jackknife enabled: ', cfg%njk, ' regions'
  end subroutine allocate_result_arrays

  subroutine deallocate_arrays()
    implicit none
    if (allocated(N2)) deallocate(N2)
    if (allocated(N3)) deallocate(N3)
    if (allocated(N2jk)) deallocate(N2jk)
    if (allocated(N3jk)) deallocate(N3jk)
    if (allocated(region)) deallocate(region)
    if (allocated(N4)) deallocate(N4)
    if (allocated(R4)) deallocate(R4)
    if (allocated(bintable6)) deallocate(bintable6)
    if (allocated(canon_bins_4pcf)) deallocate(canon_bins_4pcf)
    if (allocated(orbit_mult_4pcf)) deallocate(orbit_mult_4pcf)
    if (allocated(bintable)) deallocate(bintable)
    if (allocated(zeta3_internal)) deallocate(zeta3_internal)
    if (allocated(N4jk)) deallocate(N4jk)
    if (allocated(R4jk)) deallocate(R4jk)
    if (allocated(DD_2pcf)) deallocate(DD_2pcf)
    if (allocated(RR_2pcf)) deallocate(RR_2pcf)
    if (allocated(DD_2pcf_jk)) deallocate(DD_2pcf_jk)
    if (allocated(RR_2pcf_jk)) deallocate(RR_2pcf_jk)
    if (allocated(xi_2pcf)) deallocate(xi_2pcf)
    if (allocated(sum_pair_l2)) deallocate(sum_pair_l2)
    if (allocated(sum_pair_l4)) deallocate(sum_pair_l4)
    if (allocated(xi2_2pcf)) deallocate(xi2_2pcf)
    if (allocated(xi4_2pcf)) deallocate(xi4_2pcf)
    if (allocated(weights)) deallocate(weights)
    if (allocated(radial_bins)) deallocate(radial_bins)
  end subroutine deallocate_arrays

  ! Provenance block written at the top of every output file, so results are
  ! self-documenting: version/build, timestamp, the exact command line, and
  ! the catalogue sizes.  All lines start with '#' (numpy's default comment
  ! character), as does the column-name header that follows.
  subroutine write_provenance(unit_num)
    integer, intent(in) :: unit_num
    character(len=4000) :: cmd
    integer :: dt(8)

    call get_command(cmd)
    call date_and_time(values=dt)

    write(unit_num, '(a)') '# gramsci ' // GRAMSCI_VERSION // &
                           ' (' // trim(cfg%backend) // ' build)'
    write(unit_num, '(a,i4.4,"-",i2.2,"-",i2.2,1x,i2.2,":",i2.2,":",i2.2)') &
      '# date: ', dt(1), dt(2), dt(3), dt(5), dt(6), dt(7)
    write(unit_num, '(a)') '# command: ' // trim(cmd)
    write(unit_num, '(a,i0,a,i0)') '# n_data = ', cfg%num_data, &
                                   '   n_randoms = ', cfg%num_rand
    if (cfg%analytic) then
      write(unit_num, '(a,3(f0.4,1x))') '# analytic randoms, periodic box: ', cfg%boxsize
    else if (cfg%periodic) then
      write(unit_num, '(a,3(f0.4,1x))') '# periodic box: ', cfg%boxsize
    end if
  end subroutine write_provenance

  ! Output filename for one query mode: with a single mode this is exactly
  ! -out; with combined modes each result goes to <output_file>.<mode>.
  function mode_output_file(mode) result(fname)
    character(len=*), intent(in) :: mode
    character(len=2000) :: fname
    if (cfg%multi_mode) then
      fname = trim(cfg%output_file) // '.' // trim(mode)
    else
      fname = cfg%output_file
    end if
  end function mode_output_file

  character(len=20) function str(k)
    integer, intent(in) :: k
    write(str, *) k
    str = adjustl(str)
  end function str

  ! Allocate the 4PCF jackknife touching-sum arrays (nch = 1 for -4pcf,
  ! 2 for -4pcfp).  Unlike N2jk/N3jk these are NOT named in an OpenMP
  ! reduction: at 4PCF configuration counts a per-thread reduction copy of a
  ! (n_configs, nch, njk) array would multiply memory by the thread count,
  ! so the query kernels update them with !$OMP ATOMIC instead.  Collisions
  ! are rare because the updates scatter over n_configs*njk slots (the same
  ! reasoning as the OpenACC build's slot-strided jackknife partials).
  subroutine allocate_4pcf_jk(nch)
    integer, intent(in) :: nch
    real(kdkind) :: gb
    allocate(N4jk(cfg%n_configs_4pcf, nch, max(cfg%njk, 1)))
    allocate(R4jk(cfg%n_configs_4pcf, nch, max(cfg%njk, 1)))
    N4jk = 0.0d0
    R4jk = 0.0d0
    if (cfg%njk > 0) then
      gb = 2.0d0 * real(cfg%n_configs_4pcf, kdkind) * real(nch, kdkind) * &
           real(cfg%njk, kdkind) * 8.0d0 / 1.0d9
      if (cfg%rank == 0) &
        print '(" 4PCF jackknife accumulators: ",f8.2," GB")', gb
      if (gb > 64.0d0) then
        print *, 'ERROR: 4PCF jackknife accumulators exceed 64 GB;'
        print *, '       reduce -njk or -nbins'
        stop 1
      end if
    end if
  end subroutine allocate_4pcf_jk

  subroutine free_4pcf_jk()
    if (allocated(N4jk)) deallocate(N4jk)
    if (allocated(R4jk)) deallocate(R4jk)
  end subroutine free_4pcf_jk

  ! Delete-one jackknife mean and error of one quantity from its njk
  ! realisations:  sigma^2 = (njk-1)/njk * sum_m (x_m - xbar)^2.
  subroutine jk_mean_sigma(vals, mean, sigma)
    real(kdkind), intent(in) :: vals(:)
    real(kdkind), intent(out) :: mean, sigma
    integer :: n
    n = size(vals)
    mean = sum(vals) / real(n, kdkind)
    sigma = sqrt(real(n - 1, kdkind) / real(n, kdkind) * sum((vals - mean)**2))
  end subroutine jk_mean_sigma

end module config_module
