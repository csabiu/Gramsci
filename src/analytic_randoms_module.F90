module analytic_randoms_module
  ! ---------------------------------------------------------------------------
  ! Analytic random counts for a periodic box (-box without -ran).
  !
  ! For uniformly distributed points in a periodic volume V the expected
  ! pair / triple / quadruple counts per bin configuration are exactly
  ! computable, so no random catalogue is needed:
  !
  !   RR(bin)          = W2 * Vshell(bin) / V
  !   RRR(a<=b<=c)     = W3 * Nperm(a,b,c) * 8 pi^2 * T(a,b,c) / V^2
  !   RRRR(config)     = W4 * orbit_mult(config) * 8 pi^2 * Q(config) / V^3
  !
  ! where the W's are the weight normalizations (elementary symmetric sums of
  ! the normalized data weights, matching the "count each unordered tuple
  ! once" convention of the graph queries),
  !
  !   T(a,b,c) = int_{s1 in a, s2 in b, s3 in c} s1 s2 s3 [triangle ineq.] ds^3
  !
  ! is the classical triangle kernel (from d3r12 d3r13 = 8 pi^2 r12 r13 r23
  ! dr12 dr13 dr23), and Q is the tetrahedron kernel
  !
  !   Q = int ds12 ds13 ds23  s12 s13 s23 [tri(123)] * G(s12,s13,s23)
  !   G = volume of { r4 : |r4-p1| in b14, |r4-p2| in b24, |r4-p3| in b34 }
  !
  ! G is evaluated in cylindrical coordinates (t, w=rho^2, phi) about the 1-2
  ! axis: the phi window is analytic (arccos), the w window and its arccos
  ! kink locations are closed-form, and t is integrated with composite
  ! Gauss-Legendre split at the closed-form annulus breakpoints.  T is exact
  ! up to Gauss-Legendre order (the integrand is piecewise polynomial and all
  ! breakpoints are split explicitly).  Q uses an escalating-order rule that
  ! stops when two successive refinements agree.
  !
  ! All formulas were validated against brute-force Monte Carlo counting of
  ! uniform points in a periodic box.
  ! ---------------------------------------------------------------------------
  use kdtree2_precision_module
  use config_module
  implicit none
  private
  public :: analytic_rr, analytic_rrr, analytic_rrr_equilateral, analytic_rrrr

  real(kdkind), parameter :: PI = 3.141592653589793238462643383279502884197d0
  real(kdkind), parameter :: EIGHT_PI2 = 8.0d0 * PI * PI

  ! Precomputed Gauss-Legendre rules on [-1,1]: node k of order n is
  ! gl_x(k,n).  Orders 1..GL_MAX are filled once by init_gl (not thread-safe;
  ! called before any OpenMP region).
  integer, parameter :: GL_MAX = 64
  real(kdkind) :: gl_x(GL_MAX, GL_MAX), gl_w(GL_MAX, GL_MAX)
  logical :: gl_ready = .false.

contains

  ! ===========================================================================
  ! Weight normalizations: elementary symmetric sums e2, e3, e4 of the
  ! normalized data weights (sum w = 1), expressed via power sums.
  ! ===========================================================================
  function weight_pair() result(w)
    real(kdkind) :: w
    w = (1.0d0 - cfg%sw2) / 2.0d0
  end function weight_pair

  function weight_triple() result(w)
    real(kdkind) :: w
    w = (1.0d0 - 3.0d0*cfg%sw2 + 2.0d0*cfg%sw3) / 6.0d0
  end function weight_triple

  function weight_quad() result(w)
    real(kdkind) :: w
    w = (1.0d0 - 6.0d0*cfg%sw2 + 3.0d0*cfg%sw2*cfg%sw2 &
         + 8.0d0*cfg%sw3 - 6.0d0*cfg%sw4) / 24.0d0
  end function weight_quad

  function box_volume() result(v)
    real(kdkind) :: v
    v = cfg%boxsize(1) * cfg%boxsize(2) * cfg%boxsize(3)
  end function box_volume

  ! ===========================================================================
  ! RR: shell volume fractions.
  ! ===========================================================================
  subroutine analytic_rr(rr)
    real(kdkind), intent(out) :: rr(cfg%nbins)
    integer :: b
    real(kdkind) :: w2, v
    w2 = weight_pair()
    v = box_volume()
    do b = 1, cfg%nbins
      rr(b) = w2 * (4.0d0*PI/3.0d0) * (radial_bins(b+1)**3 - radial_bins(b)**3) / v
    end do
  end subroutine analytic_rr

  ! ===========================================================================
  ! RRR for all 3PCF configs, indexed like N3 (via bintable).
  ! ===========================================================================
  subroutine analytic_rrr(rrr)
    real(kdkind), intent(out) :: rrr(cfg%config_bins)
    integer :: i, j, k, nperm
    real(kdkind) :: w3, v2, t
    call init_gl()
    w3 = weight_triple()
    v2 = box_volume()**2
    rrr = 0.0d0
    !$OMP PARALLEL DO schedule(dynamic) private(i, j, k, nperm, t)
    do i = 1, cfg%nbins
      do j = i, cfg%nbins
        do k = j, cfg%nbins
          if (i == j .and. j == k) then
            nperm = 1
          else if (i == j .or. j == k) then
            nperm = 3
          else
            nperm = 6
          end if
          t = t_kernel(radial_bins(i), radial_bins(i+1), &
                       radial_bins(j), radial_bins(j+1), &
                       radial_bins(k), radial_bins(k+1))
          rrr(bintable(i, j, k, 1)) = w3 * nperm * EIGHT_PI2 * t / v2
        end do
      end do
    end do
    !$OMP END PARALLEL DO
  end subroutine analytic_rrr

  ! RRR for the equilateral-only mode: config (l,l,l) per radial bin.
  subroutine analytic_rrr_equilateral(rrr)
    real(kdkind), intent(out) :: rrr(cfg%nbins)
    integer :: l
    real(kdkind) :: w3, v2
    call init_gl()
    w3 = weight_triple()
    v2 = box_volume()**2
    !$OMP PARALLEL DO schedule(dynamic) private(l)
    do l = 1, cfg%nbins
      rrr(l) = w3 * EIGHT_PI2 * t_kernel(radial_bins(l), radial_bins(l+1), &
                                         radial_bins(l), radial_bins(l+1), &
                                         radial_bins(l), radial_bins(l+1)) / v2
    end do
    !$OMP END PARALLEL DO
  end subroutine analytic_rrr_equilateral

  ! ===========================================================================
  ! RRRR for all 4PCF configs (canonical order of canon_bins_4pcf).
  ! ===========================================================================
  subroutine analytic_rrrr(rrrr)
    real(kdkind), intent(out) :: rrrr(cfg%n_configs_4pcf)
    integer :: c, lvl
    integer :: bins(6)
    real(kdkind) :: w4, v3, q_prev, q
    real(kdkind) :: start, finish
    call init_gl()
    w4 = weight_quad()
    v3 = box_volume()**3
    if (cfg%rank == 0) then
      print '("Computing analytic RRRR for ",i8," configurations...")', cfg%n_configs_4pcf
    end if
    call cpu_time(start)
    !$OMP PARALLEL DO schedule(dynamic) private(c, lvl, bins, q_prev, q)
    do c = 1, cfg%n_configs_4pcf
      bins = canon_bins_4pcf(:, c)
      ! Escalate the quadrature order until two successive levels agree.
      q_prev = q_kernel(bins, 1)
      q = q_prev
      do lvl = 2, 3
        q = q_kernel(bins, lvl)
        if (abs(q - q_prev) <= 1.0d-4 * max(abs(q), 1.0d-300)) exit
        q_prev = q
      end do
      rrrr(c) = w4 * orbit_mult_4pcf(c) * EIGHT_PI2 * q / v3
    end do
    !$OMP END PARALLEL DO
    call cpu_time(finish)
    if (cfg%rank == 0) then
      print '("Analytic RRRR done (",f10.2," cpu-seconds).")', finish - start
    end if
  end subroutine analytic_rrrr

  ! ===========================================================================
  ! Gauss-Legendre tables (Newton iteration on Legendre polynomials).
  ! ===========================================================================
  subroutine init_gl()
    integer :: n, k, it
    real(kdkind) :: x, p, dp, dx
    if (gl_ready) return
    gl_x(1, 1) = 0.0d0
    gl_w(1, 1) = 2.0d0
    do n = 2, GL_MAX
      do k = 1, (n + 1) / 2
        x = cos(PI * (k - 0.25d0) / (n + 0.5d0))
        do it = 1, 100
          call legendre_pdp(n, x, p, dp)
          dx = p / dp
          x = x - dx
          if (abs(dx) < 1.0d-15) exit
        end do
        call legendre_pdp(n, x, p, dp)
        gl_x(k, n) = -x
        gl_x(n + 1 - k, n) = x
        gl_w(k, n) = 2.0d0 / ((1.0d0 - x * x) * dp * dp)
        gl_w(n + 1 - k, n) = gl_w(k, n)
      end do
    end do
    gl_ready = .true.
  end subroutine init_gl

  ! Legendre polynomial P_n(x) and derivative via the standard recurrence.
  subroutine legendre_pdp(n, x, p, dp)
    integer, intent(in) :: n
    real(kdkind), intent(in) :: x
    real(kdkind), intent(out) :: p, dp
    integer :: m
    real(kdkind) :: p0, p1
    p0 = 1.0d0
    p1 = x
    do m = 2, n
      p = ((2*m - 1) * x * p1 - (m - 1) * p0) / m
      p0 = p1
      p1 = p
    end do
    p = p1
    dp = n * (x * p1 - p0) / (x * x - 1.0d0)
  end subroutine legendre_pdp

  ! Sorted unique interior breakpoints of (lo, hi) from cand(1:nc), plus the
  ! endpoints.  pts(1:npts) with npts >= 2.
  subroutine split_points(lo, hi, cand, nc, pts, npts)
    real(kdkind), intent(in) :: lo, hi, cand(:)
    integer, intent(in) :: nc
    real(kdkind), intent(out) :: pts(:)
    integer, intent(out) :: npts
    integer :: i, j
    real(kdkind) :: eps, c
    eps = 1.0d-12 * (hi - lo)
    npts = 2
    pts(1) = lo
    pts(2) = hi
    do i = 1, nc
      c = cand(i)
      if (c <= lo + eps .or. c >= hi - eps) cycle
      ! insertion into sorted list, skipping near-duplicates
      do j = 1, npts
        if (abs(c - pts(j)) <= eps) exit
        if (c < pts(j)) then
          pts(j+1:npts+1) = pts(j:npts)
          pts(j) = c
          npts = npts + 1
          exit
        end if
      end do
    end do
  end subroutine split_points

  ! ===========================================================================
  ! T(a,b,c): triangle kernel with exact breakpoint splitting.
  ! ===========================================================================
  function t_kernel(alo, ahi, blo, bhi, clo, chi) result(tot)
    real(kdkind), intent(in) :: alo, ahi, blo, bhi, clo, chi
    real(kdkind) :: tot
    real(kdkind) :: cand(14), pts(20)
    integer :: npts, p, g
    real(kdkind) :: mid, half, x, w
    integer, parameter :: NG1 = 12
    cand = [clo - blo, clo - bhi, clo + blo, clo + bhi, &
            chi - blo, chi - bhi, chi + blo, chi + bhi, &
            blo - clo, bhi - clo, blo - chi, bhi - chi, blo, bhi]
    tot = 0.0d0
    call split_points(alo, ahi, cand, 14, pts, npts)
    do p = 1, npts - 1
      mid = 0.5d0 * (pts(p) + pts(p+1))
      half = 0.5d0 * (pts(p+1) - pts(p))
      do g = 1, NG1
        x = mid + half * gl_x(g, NG1)
        w = half * gl_w(g, NG1)
        tot = tot + w * x * t_inner(x, blo, bhi, clo, chi)
      end do
    end do
  end function t_kernel

  ! H(s1) = int_{s2 in b} s2 * K(s1,s2) ds2 with kinks split exactly.
  function t_inner(s1, blo, bhi, clo, chi) result(tot)
    real(kdkind), intent(in) :: s1, blo, bhi, clo, chi
    real(kdkind) :: tot
    real(kdkind) :: cand(7), pts(12)
    integer :: npts, p, g
    real(kdkind) :: mid, half, x, w, m1, m2
    integer, parameter :: NG2 = 8
    cand = [chi - s1, s1 - clo, s1 + clo, s1 - chi, s1 + chi, clo - s1, s1]
    tot = 0.0d0
    call split_points(blo, bhi, cand, 7, pts, npts)
    do p = 1, npts - 1
      mid = 0.5d0 * (pts(p) + pts(p+1))
      half = 0.5d0 * (pts(p+1) - pts(p))
      do g = 1, NG2
        x = mid + half * gl_x(g, NG2)
        w = half * gl_w(g, NG2)
        m1 = max(clo, abs(s1 - x))
        m2 = min(chi, s1 + x)
        if (m2 > m1) tot = tot + w * x * 0.5d0 * (m2 * m2 - m1 * m1)
      end do
    end do
  end function t_inner

  ! ===========================================================================
  ! Q(config): 5D tetrahedron kernel.
  ! bins = canonical (b12, b13, b14, b23, b24, b34); level = quadrature level.
  ! ===========================================================================
  function q_kernel(bins, level) result(tot)
    integer, intent(in) :: bins(6), level
    real(kdkind) :: tot
    integer :: n12, n13, n23, ntp, ntg, nwg
    integer :: g12, g13, g23
    real(kdkind) :: b12lo, b12hi, b13lo, b13hi, b14lo, b14hi
    real(kdkind) :: b23lo, b23hi, b24lo, b24hi, b34lo, b34hi
    real(kdkind) :: m12, h12, m13, h13, m23, h23
    real(kdkind) :: x12, w12, x13, w13, x23, w23, lo23, hi23, g

    select case (level)
    case (1)
      n12 = 8;  n13 = 8;  n23 = 10; ntp = 16; ntg = 8; nwg = 8
    case (2)
      n12 = 12; n13 = 12; n23 = 16; ntp = 32; ntg = 8; nwg = 8
    case default
      n12 = 16; n13 = 16; n23 = 20; ntp = 64; ntg = 10; nwg = 10
    end select

    b12lo = radial_bins(bins(1)); b12hi = radial_bins(bins(1) + 1)
    b13lo = radial_bins(bins(2)); b13hi = radial_bins(bins(2) + 1)
    b14lo = radial_bins(bins(3)); b14hi = radial_bins(bins(3) + 1)
    b23lo = radial_bins(bins(4)); b23hi = radial_bins(bins(4) + 1)
    b24lo = radial_bins(bins(5)); b24hi = radial_bins(bins(5) + 1)
    b34lo = radial_bins(bins(6)); b34hi = radial_bins(bins(6) + 1)

    tot = 0.0d0
    m12 = 0.5d0 * (b12lo + b12hi); h12 = 0.5d0 * (b12hi - b12lo)
    m13 = 0.5d0 * (b13lo + b13hi); h13 = 0.5d0 * (b13hi - b13lo)
    do g12 = 1, n12
      x12 = m12 + h12 * gl_x(g12, n12)
      w12 = h12 * gl_w(g12, n12)
      do g13 = 1, n13
        x13 = m13 + h13 * gl_x(g13, n13)
        w13 = h13 * gl_w(g13, n13)
        ! clip s23 to the triangle-inequality window
        lo23 = max(b23lo, abs(x12 - x13))
        hi23 = min(b23hi, x12 + x13)
        if (hi23 <= lo23) cycle
        m23 = 0.5d0 * (lo23 + hi23); h23 = 0.5d0 * (hi23 - lo23)
        do g23 = 1, n23
          x23 = m23 + h23 * gl_x(g23, n23)
          w23 = h23 * gl_w(g23, n23)
          g = g_volume(x12, x13, x23, b14lo, b14hi, b24lo, b24hi, &
                       b34lo, b34hi, ntp, ntg, nwg)
          tot = tot + w12 * w13 * w23 * x12 * x13 * x23 * g
        end do
      end do
    end do
  end function q_kernel

  ! Volume of { r4 : |r4-p1| in [p14lo,p14hi], |r4-p2| in [p24lo,p24hi],
  !                  |r4-p3| in [p34lo,p34hi] } for the triangle
  ! p1=(0,0,0), p2=(s12,0,0), p3=(x3,y3,0).
  function g_volume(s12, s13, s23, p14lo, p14hi, p24lo, p24hi, &
                    p34lo, p34hi, ntp, ntg, nwg) result(tot)
    real(kdkind), intent(in) :: s12, s13, s23
    real(kdkind), intent(in) :: p14lo, p14hi, p24lo, p24hi, p34lo, p34hi
    integer, intent(in) :: ntp, ntg, nwg
    real(kdkind) :: tot
    real(kdkind) :: x3, y3sq, y3, a1, a2, bb1, bb2, c1, c2
    real(kdkind) :: t_lo, t_hi, cand(12), pts(16)
    integer :: npts, p, pan, npan, g
    real(kdkind) :: plo, phi_, mid, half, t, wt

    tot = 0.0d0
    x3 = (s12 * s12 + s13 * s13 - s23 * s23) / (2.0d0 * s12)
    y3sq = s13 * s13 - x3 * x3
    if (y3sq <= 0.0d0) return
    y3 = sqrt(y3sq)
    a1 = p14lo * p14lo; a2 = p14hi * p14hi
    bb1 = p24lo * p24lo; bb2 = p24hi * p24hi
    c1 = p34lo * p34lo; c2 = p34hi * p34hi

    t_lo = max(-p14hi, s12 - p24hi)
    t_hi = min(p14hi, s12 + p24hi)
    if (t_hi <= t_lo) return

    cand(1) = (a1 - bb1 + s12 * s12) / (2.0d0 * s12)
    cand(2) = (a2 - bb2 + s12 * s12) / (2.0d0 * s12)
    cand(3) = (a1 - bb2 + s12 * s12) / (2.0d0 * s12)
    cand(4) = (a2 - bb1 + s12 * s12) / (2.0d0 * s12)
    cand(5) = p14lo;  cand(6) = -p14lo
    cand(7) = p14hi;  cand(8) = -p14hi
    cand(9) = s12 + p24lo;  cand(10) = s12 - p24lo
    cand(11) = s12 + p24hi; cand(12) = s12 - p24hi

    call split_points(t_lo, t_hi, cand, 12, pts, npts)
    do p = 1, npts - 1
      npan = max(1, ceiling(real(ntp, kdkind) * (pts(p+1) - pts(p)) / (t_hi - t_lo)))
      do pan = 1, npan
        plo = pts(p) + (pts(p+1) - pts(p)) * (pan - 1) / real(npan, kdkind)
        phi_ = pts(p) + (pts(p+1) - pts(p)) * pan / real(npan, kdkind)
        mid = 0.5d0 * (plo + phi_)
        half = 0.5d0 * (phi_ - plo)
        do g = 1, ntg
          t = mid + half * gl_x(g, ntg)
          wt = half * gl_w(g, ntg)
          tot = tot + wt * g_w_integral(t, s12, x3, y3, y3sq, a1, a2, bb1, bb2, c1, c2, nwg)
        end do
      end do
    end do
  end function g_volume

  ! Inner integral over w = rho^2 at fixed t:  int dw/2 * Phi(t, w).
  function g_w_integral(t, s12, x3, y3, y3sq, a1, a2, bb1, bb2, c1, c2, nwg) result(tot)
    real(kdkind), intent(in) :: t, s12, x3, y3, y3sq, a1, a2, bb1, bb2, c1, c2
    integer, intent(in) :: nwg
    real(kdkind) :: tot
    real(kdkind) :: wlo, whi, e, cand(6), pts(10)
    real(kdkind) :: disc, r, gg, cc
    integer :: nc, npts, p, g, icase
    real(kdkind) :: mid, half, w, wt, d, rt, z1, z2

    tot = 0.0d0
    wlo = max(a1 - t * t, bb1 - (t - s12)**2, 0.0d0)
    whi = min(a2 - t * t, bb2 - (t - s12)**2)
    if (whi <= wlo) return

    e = (t - x3)**2 + y3sq
    nc = 0
    do icase = 1, 2
      if (icase == 1) then
        cc = c1
      else
        cc = c2
      end if
      disc = y3sq - e + cc
      if (disc > 0.0d0) then
        r = sqrt(disc)
        gg = y3 + r
        if (gg > 0.0d0) then
          nc = nc + 1; cand(nc) = gg * gg
        end if
        gg = y3 - r
        if (gg > 0.0d0) then
          nc = nc + 1; cand(nc) = gg * gg
        end if
        gg = -y3 + r
        if (gg > 0.0d0) then
          nc = nc + 1; cand(nc) = gg * gg
        end if
      end if
    end do

    call split_points(wlo, whi, cand, nc, pts, npts)
    do p = 1, npts - 1
      mid = 0.5d0 * (pts(p) + pts(p+1))
      half = 0.5d0 * (pts(p+1) - pts(p))
      do g = 1, nwg
        w = mid + half * gl_x(g, nwg)
        wt = half * gl_w(g, nwg)
        if (w <= 0.0d0) cycle
        d = (t - x3)**2 + w + y3sq
        rt = 2.0d0 * y3 * sqrt(w)
        z1 = max(-1.0d0, min(1.0d0, (d - c2) / rt))
        z2 = max(-1.0d0, min(1.0d0, (d - c1) / rt))
        tot = tot + wt * 0.5d0 * 2.0d0 * (acos(z1) - acos(z2))
      end do
    end do
  end function g_w_integral

end module analytic_randoms_module
