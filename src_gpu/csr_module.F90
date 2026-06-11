! ---------------------------------------------------------------------------
! CSR (Compressed Sparse Row) adjacency representation for GPU offloading.
!
! The original graph stores adjacency lists as an array of derived types
! (output(:)%id, output(:)%dist, etc.) with allocatable per-node arrays.
! That jagged structure cannot be mapped to GPU memory.  This module flattens
! it into four contiguous arrays:
!
!   csr_ptr(i) .. csr_ptr(i+1)-1   — range of edge indices for node i
!   csr_id(e)                       — neighbor node index for edge e
!   csr_dist(e)                     — distance-bin index (int8)
!   csr_mu(e)                       — mu-bin index (int8)
!   csr_phi(e)                       — direction-pixel index (int8, 4PCFp only)
!
! call build_csr()           after create_graph()
! call deallocate_csr()      at cleanup
! ---------------------------------------------------------------------------
module csr_module
  use iso_fortran_env, only: int8, int64
  use iso_c_binding, only: c_int, c_size_t
  use config_module
  implicit none

  ! glibc malloc_trim: releases freed heap pages back to the OS (madvise).
  ! Needed in build_csr — without it the freed jagged-graph rows stay in
  ! glibc's free lists while the CSR arrays commit, and RSS peaks at
  ! jagged+CSR (≈93 GB for DESI LRG at rmax=150 → OOM on a 125 GB box).
  interface
    function malloc_trim(pad) bind(c, name='malloc_trim')
      import :: c_int, c_size_t
      integer(c_size_t), value :: pad
      integer(c_int) :: malloc_trim
    end function malloc_trim
  end interface

  ! csr_ptr / csr_total_edges are 64-bit: large catalogs (e.g. DESI LRG at
  ! rmax>~90) exceed 2^31-1 edges.  Node ids and per-node counts stay int32.
  integer(int64), allocatable :: csr_ptr(:) ! shape (N+1), 1-based offsets
  integer, allocatable :: csr_id(:)         ! shape (total_edges)
  integer(int8), allocatable :: csr_dist(:) ! shape (total_edges)
  integer(int8), allocatable :: csr_mu(:)   ! shape (total_edges)
  integer(int8), allocatable :: csr_phi(:)  ! shape (total_edges), 4PCFp only
  integer(int64) :: csr_total_edges = 0

  ! Device-memory sizing for the chunked GPU kernels.
  ! BYTES_PER_EDGE: csr_id (int32) + csr_dist (int8) per resident edge window.
  ! RESERVE_BYTES: csr_ptr + weights + buffer + partial accumulators + CUDA
  ! context — generously over-budgeted.
  integer(int64), parameter :: BYTES_PER_EDGE = 5_int64
  integer(int64), parameter :: RESERVE_BYTES  = 1500000000_int64

contains

  ! Build CSR arrays from the node-list representation in output(:).
  ! Call this immediately after create_graph().  Each jagged row is FREED as
  ! soon as it is copied, so peak host RAM is ~one full copy of the graph,
  ! not two (matters for billion-edge runs).  output(:) itself is deallocated
  ! here as well — do not touch it after this call.
  subroutine build_csr()
    integer :: i, k, N, rc
    integer(int64) :: base
    logical :: want_mu

    N = cfg%num_data + cfg%num_rand

    ! mu bins are only consumed by (CPU) RSD queries; skip the array otherwise
    ! — it is 1 byte/edge of dead weight for isotropic runs.
    want_mu = (cfg%nmu > 1)

    ! Pass 1: prefix sums → csr_ptr (64-bit; can exceed 2^31-1 edges)
    allocate(csr_ptr(N + 1))
    csr_ptr(1) = 1
    do i = 1, N
      csr_ptr(i + 1) = csr_ptr(i) + output(i)%nn
    end do
    csr_total_edges = csr_ptr(N + 1) - 1

    ! Allocate flat edge arrays
    allocate(csr_id(csr_total_edges))
    allocate(csr_dist(csr_total_edges))
    if (want_mu) allocate(csr_mu(csr_total_edges))
    if (cfg%four_pcf_parity) allocate(csr_phi(csr_total_edges))

    ! Pass 2: copy from jagged output(:), freeing each row as we go
    do i = 1, N
      base = csr_ptr(i) - 1
      do k = 1, output(i)%nn
        csr_id  (base + k) = output(i)%id  (k)
        csr_dist(base + k) = output(i)%dist(k)
      end do
      if (want_mu) then
        do k = 1, output(i)%nn
          csr_mu(base + k) = output(i)%mu(k)
        end do
      end if
      if (cfg%four_pcf_parity .and. allocated(output(i)%phi)) then
        do k = 1, output(i)%nn
          csr_phi(base + k) = output(i)%phi(k)
        end do
      end if
      call output(i)%destroy()
      ! Hand freed row pages back to the OS so RSS tracks the live data,
      ! not the high-water mark (see malloc_trim interface note above).
      if (mod(i, 131072) == 0) rc = malloc_trim(0_c_size_t)
    end do
    deallocate(output)
    rc = malloc_trim(0_c_size_t)

    if (cfg%rank == 0) &
      print '("CSR graph: ",i0," nodes, ",i0," edges")', N, csr_total_edges
  end subroutine build_csr

  ! ---------------------------------------------------------------------------
  ! Binary search for the distance bin of edge (from_node → to_node) in the
  ! CSR half-graph.  Returns 0_int8 if the edge is absent.
  !
  ! !$ACC ROUTINE SEQ makes this callable from inside GPU parallel regions.
  ! The arrays csr_ptr / csr_id / csr_dist must be in the caller's ACC DATA
  ! region so the device copy is used here.
  ! ---------------------------------------------------------------------------
  function csr_find_dist_bin(ptr, id_arr, dist_arr, from_node, to_node) result(bin)
    !$ACC ROUTINE SEQ
    integer(int64), intent(in) :: ptr(:)
    integer,        intent(in) :: id_arr(:)
    integer(int8),  intent(in) :: dist_arr(:)
    ! VALUE: by-reference scalars passed into device routines get a one-shot
    ! present-table copy that is NOT refreshed on later kernel launches —
    ! a stale-argument hazard when the caller updates them between launches.
    integer, value :: from_node, to_node
    integer(int8) :: bin
    integer :: lo, hi, mid, nn
    integer(int64) :: base

    bin = 0_int8
    base = ptr(from_node) - 1
    nn   = int(ptr(from_node + 1) - ptr(from_node))
    if (nn == 0) return

    ! Early exit: to_node outside the sorted range
    if (id_arr(base + 1) > to_node .or. id_arr(base + nn) < to_node) return

    lo = 1
    hi = nn
    do while (lo <= hi)
      mid = lo + (hi - lo) / 2
      if (id_arr(base + mid) == to_node) then
        bin = dist_arr(base + mid)
        return
      else if (id_arr(base + mid) < to_node) then
        lo = mid + 1
      else
        hi = mid - 1
      end if
    end do
  end function csr_find_dist_bin

  ! ---------------------------------------------------------------------------
  ! Windowed variant of csr_find_dist_bin for chunked GPU processing: the edge
  ! arrays passed in hold only the window [off+1 : off+size(id_win)] of the
  ! full CSR, so global edge index e maps to local index e - off.  The caller
  ! guarantees from_node's row lies entirely inside the window.
  ! ---------------------------------------------------------------------------
  function csr_find_dist_bin_win(ptr, id_win, dist_win, off, from_node, to_node) result(bin)
    !$ACC ROUTINE SEQ
    integer(int64), intent(in) :: ptr(:)
    integer,        intent(in) :: id_win(:)
    integer(int8),  intent(in) :: dist_win(:)
    ! VALUE: see csr_find_dist_bin — by-reference scalars go stale across
    ! kernel launches when the host updates them (off changes per tile).
    integer(int64), value :: off
    integer,        value :: from_node, to_node
    integer(int8) :: bin
    integer :: lo, hi, mid, nn
    integer(int64) :: base

    bin = 0_int8
    base = ptr(from_node) - 1 - off
    nn   = int(ptr(from_node + 1) - ptr(from_node))
    if (nn == 0) return

    if (id_win(base + 1) > to_node .or. id_win(base + nn) < to_node) return

    lo = 1
    hi = nn
    do while (lo <= hi)
      mid = lo + (hi - lo) / 2
      if (id_win(base + mid) == to_node) then
        bin = dist_win(base + mid)
        return
      else if (id_win(base + mid) < to_node) then
        lo = mid + 1
      else
        hi = mid - 1
      end if
    end do
  end function csr_find_dist_bin_win

  ! ---------------------------------------------------------------------------
  ! Runtime GPU sizing helpers
  ! ---------------------------------------------------------------------------

  ! Free device memory in bytes; -1 when running without a GPU (host fallback).
  function gpu_free_mem_bytes() result(bytes)
    use openacc
    integer(int64) :: bytes
    integer(acc_device_kind) :: devtype
    integer :: devnum

    devtype = acc_get_device_type()
    if (devtype == acc_device_host .or. devtype == acc_device_none) then
      bytes = -1_int64
      return
    end if
    devnum = acc_get_device_num(devtype)
    bytes = int(acc_get_property(devnum, devtype, acc_property_free_memory), int64)
    if (bytes <= 0) bytes = -1_int64
  end function gpu_free_mem_bytes

  ! Longest adjacency row in the CSR graph (sizes per-hub scratch like lmat).
  function csr_max_row_len() result(mx)
    integer :: mx, i, n
    mx = 0
    n = cfg%num_data + cfg%num_rand
    do i = 1, n
      mx = max(mx, int(csr_ptr(i + 1) - csr_ptr(i)))
    end do
  end function csr_max_row_len

  ! Edge capacity of one device window when `nresident` windows plus
  ! `extra_bytes` of kernel scratch must fit alongside the fixed overhead.
  ! Override with GRAMSCI_GPU_WIN_EDGES (testing / unusual devices).
  function csr_edge_window(nresident, extra_bytes) result(wedges)
    integer, intent(in) :: nresident
    integer(int64), intent(in) :: extra_bytes
    integer(int64) :: wedges, freeb, avail
    character(32) :: env
    integer :: stat

    call get_environment_variable('GRAMSCI_GPU_WIN_EDGES', env, status=stat)
    if (stat == 0) then
      read(env, *) wedges
      return
    end if

    freeb = gpu_free_mem_bytes()
    if (freeb < 0) then
      ! Host fallback: no device limit — single window.
      wedges = max(csr_total_edges, 1_int64)
      return
    end if

    avail = freeb - RESERVE_BYTES - extra_bytes
    wedges = avail / int(nresident, int64) / BYTES_PER_EDGE
    if (wedges < 50000000_int64) then
      print '("ERROR: insufficient device memory: ",i0," B free, ",i0," B scratch")', &
            freeb, extra_bytes
      stop
    end if
  end function csr_edge_window

  ! Row-aligned window boundaries: window w covers node rows
  ! [split_id(w), split_id(w+1)-1] = edges [split_e(w), split_e(w+1)-1],
  ! each window holding at most win_edges edges.
  subroutine csr_make_splits(win_edges, nwin, split_id, split_e)
    integer(int64), intent(in) :: win_edges
    integer, intent(out) :: nwin
    integer, allocatable, intent(out) :: split_id(:)
    integer(int64), allocatable, intent(out) :: split_e(:)
    integer :: n_nodes, w, lo, hi, midp
    integer(int64) :: tgt

    n_nodes = cfg%num_data + cfg%num_rand
    nwin = int((csr_total_edges - 1) / win_edges) + 1

    allocate(split_id(nwin + 1), split_e(nwin + 1))
    split_id(1) = 1
    split_e(1)  = 1
    do w = 1, nwin
      if (w == nwin) then
        split_id(w + 1) = n_nodes + 1
        split_e(w + 1)  = csr_total_edges + 1
      else
        ! Largest node id whose row still starts within the edge budget.
        tgt = split_e(w) + win_edges
        lo = split_id(w)
        hi = n_nodes + 1
        do while (lo < hi)
          midp = (lo + hi + 1) / 2
          if (csr_ptr(midp) <= tgt) then
            lo = midp
          else
            hi = midp - 1
          end if
        end do
        split_id(w + 1) = lo
        split_e(w + 1)  = csr_ptr(lo)
      end if
    end do
  end subroutine csr_make_splits

  subroutine deallocate_csr()
    if (allocated(csr_ptr))  deallocate(csr_ptr)
    if (allocated(csr_id))   deallocate(csr_id)
    if (allocated(csr_dist)) deallocate(csr_dist)
    if (allocated(csr_mu))   deallocate(csr_mu)
    if (allocated(csr_phi))  deallocate(csr_phi)
    csr_total_edges = 0
  end subroutine deallocate_csr

end module csr_module
