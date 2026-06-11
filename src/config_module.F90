module config_module
  use kdtree2_precision_module
  use node_module
  use iso_fortran_env, only: int8
  implicit none

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
    logical :: cut = .false.
    logical :: DOMPI = .false.
    logical :: rancat = .false.
    logical :: two_pcf = .false.
    logical :: three_pcf = .false.
    logical :: three_pcf_eq = .false.
    logical :: four_pcf = .false.
    logical :: four_pcf_parity = .false.
    logical :: half_graph = .true.
    integer :: n_configs_4pcf = 0
    ! Direction pixel parameters for 4PCF parity (n_theta * n_phi <= 127 for int8)
    integer :: n_theta_dir = 4
    integer :: n_phi_dir = 16
    integer :: n_dir_pixels = 64   ! = n_theta_dir * n_phi_dir
    logical :: loadran = .false.
    logical :: saveran = .false.
    character(len=2000) :: file1 = '', file2 = '', ranfile = ''
    character(len=2000) :: output_file = 'result.txt'
    ! MPI state
    integer :: rank = 0
    integer :: num_tasks = 1
    integer :: master = 0
    ! Data dimensions
    integer :: num_data = 0
    integer :: num_rand = 0
  end type gramsci_config

  ! Singleton config instance
  type(gramsci_config), save :: cfg

  ! Shared data arrays (accessible to all modules)
  type(node), dimension(:), allocatable :: output
  real(kdkind), dimension(:,:), allocatable :: points
  real(kdkind), allocatable :: weights(:), radial_bins(:)
  integer, allocatable :: buffer(:)
  real(kdkind), allocatable :: N2(:,:,:), N3(:,:,:)
  integer, allocatable :: bintable(:,:,:,:)
  ! 4PCF parity arrays
  real(kdkind), allocatable :: N4(:,:), R4(:,:)
  integer, allocatable :: bintable6(:,:,:,:,:,:)
  ! Canonical bin 6-tuple for each 4PCF config index (shape [6, n_configs_4pcf])
  integer, allocatable :: canon_bins_4pcf(:,:)
  ! Internal 2PCF for disconnected 4PCF subtraction
  real(kdkind), allocatable :: DD_2pcf(:), RR_2pcf(:), xi_2pcf(:)

contains

  subroutine default_params()
    cfg%d = 3
    cfg%cut = .false.
    cfg%logbins = .false.
    cfg%nbins = 0
    cfg%rmin = 0.0d0
    cfg%rmax = 0.0d0
    cfg%output_file = 'result.txt'
    cfg%RSD = .false.
    cfg%nmu = 1
    cfg%loadran = .false.
    cfg%saveran = .false.
    cfg%rancat = .false.
    cfg%three_pcf_eq = .false.
    cfg%three_pcf = .false.
    cfg%four_pcf = .false.
    cfg%four_pcf_parity = .false.
    cfg%two_pcf = .false.
    cfg%half_graph = .true.
  end subroutine default_params

  subroutine parseOptions()
    use extension, only: getArgument, nArguments
    implicit none
    integer(kind=4) :: i, n
    character(len=6) :: opt
    character(len=2000) :: arg

    n = nArguments()
    if (n < 6 .and. cfg%rank == cfg%master) then
      print *, ' '
      print *, 'Not enough input parameters. Please read the following help info'
      print *, ' '
      call print_help()
      stop
    end if

    i = 1
    do while (i <= n)
      call getArgument(i, opt)
      select case (opt)
      case ('-gal')
        call getArgument(i+1, arg)
        cfg%file1 = trim(arg)
        i = i + 2
      case ('-ran')
        cfg%rancat = .true.
        call getArgument(i+1, arg)
        cfg%file2 = trim(arg)
        i = i + 2
      case ('-out')
        call getArgument(i+1, arg)
        cfg%output_file = trim(arg)
        i = i + 2
      case ('-rmin')
        call getArgument(i+1, arg)
        read(arg, *) cfg%rmin
        i = i + 2
      case ('-rmax')
        call getArgument(i+1, arg)
        read(arg, *) cfg%rmax
        i = i + 2
      case ('-nbins')
        call getArgument(i+1, arg)
        read(arg, *) cfg%nbins
        i = i + 2
      case ('-nmu')
        call getArgument(i+1, arg)
        read(arg, *) cfg%nmu
        if (cfg%nmu >= 2) then
          cfg%RSD = .true.
          if (cfg%rank == cfg%master) print *, 'Anisotropic analysis requested'
        end if
        i = i + 2
      case ('-cut')
        call getArgument(i+1, arg)
        cfg%cut = .true.
        cfg%file1 = trim(arg)
        if (cfg%rank == cfg%master) print *, 'Treating input catalogue as subsample'
        if (cfg%rank == cfg%master) print *, 'will ignore -gal -ran options'
        cfg%output_file = cfg%file1
        i = i + 2
      case ('-mpi')
        cfg%DOMPI = .true.
        i = i + 1
      case ('-log')
        cfg%logbins = .true.
        if (cfg%rank == cfg%master) print *, 'Using logarithmically spaced bins'
        i = i + 1
      case ('-help')
        call print_help()
        stop
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
      case default
        print '("unknown option ",a," ignored")', trim(opt)
        stop
      end select
    end do

    ! --- Post-parse validation ---
    if (cfg%nbins <= 0) then
      print *, 'ERROR: -nbins must be specified and > 0'
      stop
    end if
    if (cfg%nbins > 127) then
      print *, 'ERROR: -nbins must be <= 127 (bin indices stored as int8)'
      stop
    end if
    if (cfg%nmu > 127) then
      print *, 'ERROR: -nmu must be <= 127 (mu indices stored as int8)'
      stop
    end if
    if (cfg%n_theta_dir * cfg%n_phi_dir > 127) then
      print *, 'ERROR: n_theta_dir * n_phi_dir must be <= 127 (phi pixel stored as int8)'
      stop
    end if
  end subroutine parseOptions

  subroutine print_help()
    print *, 'PURPOSE: Code for calculating the N-PCF of a 3D point set.'
    print *, ' '
    print *, 'CALLING SEQUENCE:'
    print *, '      gramsci [-gal galaxy file] [-ran ranfile (optional)]'
    print *, '              [-rmin Rmin] [-rmax Rmax] [-nbins Nbins]'
    print *, '              [-nmu Nmu] [-out out_file]'
    print *, '              [-2pcf | -3pcf | -equi | -4pcf | -4pcfp]'
    print *, '              [-log] [-cut]'
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
    print *, '       -cut    use file format with selection cuts'
    print *, ' '
    print *, 'QUERY MODES:'
    print *, '       -2pcf   2-point correlation function'
    print *, '       -3pcf   3-point correlation function (all configs)'
    print *, '       -equi   3-point correlation function (equilateral only)'
    print *, '       -4pcf   4-point correlation function (all configs)'
    print *, '       -4pcfp  4-point correlation function (all configs + parity)'
    print *, ' '
  end subroutine print_help

  subroutine create_binlookup()
    integer :: i, j, k

    if (cfg%three_pcf) then
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
  end subroutine allocate_result_arrays

  subroutine deallocate_arrays()
    implicit none
    if (allocated(N2)) deallocate(N2)
    if (allocated(N3)) deallocate(N3)
    if (allocated(N4)) deallocate(N4)
    if (allocated(R4)) deallocate(R4)
    if (allocated(bintable6)) deallocate(bintable6)
    if (allocated(canon_bins_4pcf)) deallocate(canon_bins_4pcf)
    if (allocated(bintable)) deallocate(bintable)
    if (allocated(DD_2pcf)) deallocate(DD_2pcf)
    if (allocated(RR_2pcf)) deallocate(RR_2pcf)
    if (allocated(xi_2pcf)) deallocate(xi_2pcf)
    if (allocated(weights)) deallocate(weights)
    if (allocated(radial_bins)) deallocate(radial_bins)
  end subroutine deallocate_arrays

  character(len=20) function str(k)
    integer, intent(in) :: k
    write(str, *) k
    str = adjustl(str)
  end function str

end module config_module
