import os
import subprocess
import tempfile
import numpy as np
from pathlib import Path

_REPO_ROOT = Path(__file__).parent.parent.parent
_DEFAULT_BINARY = _REPO_ROOT / 'bin' / 'gramsci'


def _find_binary(binary):
    if binary is None:
        binary = os.environ.get('GRAMSCI_BIN')
    if binary is None:
        if _DEFAULT_BINARY.exists():
            return str(_DEFAULT_BINARY)
        return 'gramsci'  # resolved via PATH
    binary = str(binary)
    # The subprocess runs with cwd set to a temp dir, so a relative path
    # (e.g. GRAMSCI_BIN=bin/gramsci) must be anchored to the caller's cwd.
    if os.sep in binary or os.path.exists(binary):
        return os.path.abspath(binary)
    return binary


def _prepare_randoms(randoms_pos, randoms_weights, box, what):
    """Validate the randoms/box combination and default the weights to 1."""
    if randoms_pos is None:
        if box is None:
            raise ValueError(
                f"either randoms_pos or box is required: without them the "
                f"{what} denominator is identically zero and the result "
                f"would be Inf/NaN"
            )
        return None, None
    randoms_pos = np.asarray(randoms_pos, dtype=np.float64)
    if randoms_weights is None:
        randoms_weights = np.ones(len(randoms_pos))
    return randoms_pos, randoms_weights


def _write_catalog(path, positions, weights):
    data = np.column_stack([np.asarray(positions, dtype=np.float64),
                            np.asarray(weights, dtype=np.float64)])
    np.savetxt(path, data, fmt='%.10e')


def _run(args, cwd):
    result = subprocess.run(args, capture_output=True, text=True, cwd=cwd)
    if result.returncode != 0:  # gramsci error paths exit with a nonzero code
        raise RuntimeError(
            f"gramsci exited with code {result.returncode}:\n"
            f"{result.stdout[-2000:]}\n{result.stderr[-500:]}"
        )
    return result.stdout


def _load(path, stdout=''):
    if not os.path.exists(path):
        raise FileNotFoundError(
            f"gramsci did not produce output at {path}. Binary output:\n"
            f"{stdout[-2000:]}"
        )
    # '#' marks the provenance block and column names; 'r' keeps outputs
    # written by older builds (whose header wrapped across lines) readable.
    return np.loadtxt(path, comments=('#', 'r'), ndmin=2)


def _load_jk(out, stdout, njk):
    """Load <out>.jk / <out>.jkerr written by -njk.

    Returns (real, err): the realisation table with its rows reshaped to
    (n_rows_of_main_output, njk, n_cols) and the .jkerr summary array.
    """
    err = _load(out + '.jkerr', stdout)
    real = _load(out + '.jk', stdout)
    real = real.reshape(err.shape[0], njk, real.shape[1])
    return real, err


def hartlap_factor(njk, nbins):
    """Multiplicative correction for the inverse of a jackknife covariance.

    An inverse covariance estimated from njk realisations over nbins bins
    is biased high; multiply the inverse by (njk - nbins - 2)/(njk - 1)
    (Hartlap, Simon & Schneider 2007) before using it in a likelihood.
    Requires njk > nbins + 2.
    """
    if njk <= nbins + 2:
        raise ValueError(
            f"njk = {njk} must exceed nbins + 2 = {nbins + 2} for an "
            f"invertible jackknife covariance")
    return (njk - nbins - 2) / (njk - 1)


def _jk_cov(real, col):
    """Delete-one jackknife covariance from the realisation table."""
    x = real[:, :, col]                     # (nrows, njk)
    njk = x.shape[1]
    d = x - x.mean(axis=1, keepdims=True)
    return (njk - 1) / njk * (d @ d.T)


class _JKCovMixin:
    """jk_covariance() for result objects computed with njk > 0."""

    def jk_covariance(self, column=None):
        """Delete-one jackknife covariance matrix of one estimator column.

        C_ij = (njk-1)/njk sum_m (x_i^m - xbar_i)(x_j^m - xbar_j), with
        row order matching the result arrays.  column defaults to the
        primary estimator of the mode; see the class's _jk_columns for the
        available names.  Remember hartlap_factor() when inverting.
        """
        if not hasattr(self, 'jk_real'):
            raise ValueError(
                "no jackknife realisations on this result: rerun the "
                "compute_* call with njk > 0")
        cols = self._jk_columns()
        if column is None:
            column = self._jk_default
        if column not in cols:
            raise ValueError(
                f"unknown column {column!r}; available: {sorted(cols)}")
        return _jk_cov(self.jk_real, cols[column])

    def jk_inverse_covariance(self, column=None, hartlap=True):
        """Inverse jackknife covariance, Hartlap-corrected BY DEFAULT.

        The raw inverse of a noisy covariance estimate is biased high;
        this multiplies it by (njk - nbins - 2)/(njk - 1) (Hartlap, Simon
        & Schneider 2007) so it can go straight into a Gaussian
        likelihood.  Pass hartlap=False for the uncorrected inverse.
        Requires njk > nbins + 2 (more regions than bins), otherwise the
        covariance is singular and this raises.
        """
        C = self.jk_covariance(column)
        Cinv = np.linalg.inv(C)
        if hartlap:
            Cinv = Cinv * hartlap_factor(self.njk, C.shape[0])
        return Cinv


def _box_arg(box):
    """Format the -box argument: scalar L or (Lx, Ly, Lz)."""
    box = np.atleast_1d(np.asarray(box, dtype=np.float64))
    if box.size == 1:
        return f'{box[0]:.10g}'
    if box.size == 3:
        return ','.join(f'{v:.10g}' for v in box)
    raise ValueError("box must be a scalar side length or a 3-sequence (Lx, Ly, Lz)")


# ---------------------------------------------------------------------------
# Result containers
# ---------------------------------------------------------------------------

class TwoPCFResult(_JKCovMixin):
    """Result of compute_2pcf.

    Attributes
    ----------
    r : midpoint of each radial bin
    r_min, r_max : bin edges
    mu_min, mu_max : mu-bin edges (constant when nmu = 1)
    NN, RR : weighted pair counts
    xi : correlation function estimate
    nmu : number of mu bins

    With nmu > 1 the arrays have nbins*nmu entries, ordered with the nmu
    mu bins consecutive within each radial bin.
    """
    def __init__(self, data):
        # columns: r_min r_max mu_min mu_max NN RR xi
        self.r_min  = data[:, 0]
        self.r_max  = data[:, 1]
        self.mu_min = data[:, 2]
        self.mu_max = data[:, 3]
        self.NN     = data[:, 4]
        self.RR     = data[:, 5]
        self.xi     = data[:, 6]
        self.r      = 0.5 * (self.r_min + self.r_max)
        self.nmu    = len(np.unique(self.mu_min))

    _jk_default = 'xi'

    def _jk_columns(self):
        # .jk columns: rmin rmax mumin mumax ireal NN RR xi
        return {'NN': 5, 'RR': 6, 'xi': 7}

    def __repr__(self):
        return f"TwoPCFResult(nbins={len(self.r) // self.nmu}, nmu={self.nmu})"


class ThreePCFResult(_JKCovMixin):
    """Result of compute_3pcf (all triangle configurations).

    Attributes
    ----------
    r1, r2, r3 : midpoints of triangle side bins
    r{1,2,3}_min, r{1,2,3}_max : bin edges
    mu_min, mu_max : mu-bin edges (present only when nmu > 1)
    NNN, RRR : weighted triplet counts
    zeta : 3PCF estimate  zeta = NNN / RRR
    n_configs : number of rows (triangle configs, x nmu when nmu > 1)
    """
    def __init__(self, data):
        # isotropic: 9 columns; RSD (nmu > 1): 11 columns with the mu-bin
        # edges inserted at columns 6-7
        self.r1_min = data[:, 0]; self.r1_max = data[:, 1]
        self.r2_min = data[:, 2]; self.r2_max = data[:, 3]
        self.r3_min = data[:, 4]; self.r3_max = data[:, 5]
        if data.shape[1] == 11:
            self.mu_min = data[:, 6]
            self.mu_max = data[:, 7]
            self.NNN    = data[:, 8]
            self.RRR    = data[:, 9]
            self.zeta   = data[:, 10]
        elif data.shape[1] == 9:
            self.mu_min = None
            self.mu_max = None
            self.NNN    = data[:, 6]
            self.RRR    = data[:, 7]
            self.zeta   = data[:, 8]
        else:
            raise ValueError(
                f"unexpected 3PCF output with {data.shape[1]} columns "
                "(expected 9, or 11 with nmu > 1)"
            )
        self.r1     = 0.5 * (self.r1_min + self.r1_max)
        self.r2     = 0.5 * (self.r2_min + self.r2_max)
        self.r3     = 0.5 * (self.r3_min + self.r3_max)
        self.n_configs = len(self.zeta)

    _jk_default = 'zeta'

    def _jk_columns(self):
        # .jk columns: 6 edges [+ 2 mu-bin edges with RSD] + ireal NNN RRR zeta
        base = 9 if self.jk_real.shape[2] == 12 else 7
        return {'NNN': base, 'RRR': base + 1, 'zeta': base + 2}

    def __repr__(self):
        return f"ThreePCFResult(n_configs={self.n_configs})"


class FourPCFResult(_JKCovMixin):
    """Result of compute_4pcf.

    Standard (parity=False)
    -----------------------
    NNNN, RRRR : weighted quadruplet counts
    zeta       : total 4PCF
    zeta_disc  : disconnected (lower-order) contribution
    zeta_conn  : connected 4PCF  (= zeta - zeta_disc)

    Parity-decomposed (parity=True)
    --------------------------------
    NNNN, RRRR   : even-channel counts
    zeta_even    : parity-even 4PCF
    NNNN_odd, RRRR_odd : odd-channel counts
    zeta_odd     : parity-odd 4PCF
    """
    def __init__(self, data, parity):
        self.r12_min = data[:, 0];  self.r12_max = data[:, 1]
        self.r13_min = data[:, 2];  self.r13_max = data[:, 3]
        self.r14_min = data[:, 4];  self.r14_max = data[:, 5]
        self.r23_min = data[:, 6];  self.r23_max = data[:, 7]
        self.r24_min = data[:, 8];  self.r24_max = data[:, 9]
        self.r34_min = data[:, 10]; self.r34_max = data[:, 11]
        self.n_configs = len(data)
        self.parity = parity

        if parity:
            self.NNNN      = data[:, 12]
            self.RRRR      = data[:, 13]
            self.zeta_even = data[:, 14]
            self.NNNN_odd  = data[:, 15]
            self.RRRR_odd  = data[:, 16]
            self.zeta_odd  = data[:, 17]
        else:
            self.NNNN      = data[:, 12]
            self.RRRR      = data[:, 13]
            self.zeta      = data[:, 14]
            self.zeta_disc = data[:, 15]
            self.zeta_conn = data[:, 16]

    @property
    def _jk_default(self):
        return 'zeta_conn_even' if self.parity else 'zeta_conn'

    def _jk_columns(self):
        # .jk columns: 12 edges + ireal + counts/estimators (parity adds odd)
        if self.parity:
            return {'NNNN': 13, 'RRRR': 14, 'zeta_even': 15, 'NNNN_odd': 16,
                    'zeta_odd': 17, 'zeta_disc': 18, 'zeta_conn_even': 19,
                    'zeta_conn_odd': 20}
        return {'NNNN': 13, 'RRRR': 14, 'zeta': 15, 'zeta_disc': 16,
                'zeta_conn': 17}

    def __repr__(self):
        kind = 'parity' if self.parity else 'connected'
        return f"FourPCFResult(n_configs={self.n_configs}, kind={kind!r})"


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------

def compute_2pcf(positions, weights, randoms_pos=None, randoms_weights=None,
                 rmin=1.0, rmax=30.0, nbins=10, nmu=1, box=None, binary=None,
                 njk=0):
    """Compute the 2-point correlation function.

    Parameters
    ----------
    positions : array (N, 3)  — x, y, z coordinates
    weights   : array (N,)    — per-particle weights
    randoms_pos, randoms_weights : optional random catalog for RR estimation
    rmin, rmax : float        — radial range
    nbins      : int          — number of radial bins
    nmu        : int          — number of mu bins (1 = isotropic)
    box        : float or (3,) — periodic box side length(s).  Separations use
                 the minimum image; if no randoms are given, RR is computed
                 analytically and no random catalog is needed.  Requires
                 rmax < L/2 and nmu = 1.
    binary     : str or Path  — path to gramsci binary (auto-detected if None)
    njk        : int — number of delete-one jackknife regions (0 = off).
                 Regions are equal-count ANGULAR patches of the sky as seen
                 from the origin, so radial and angular systematics are never
                 mixed.  Requires randoms; incompatible with box.  Adds
                 .xi_jk_mean, .xi_jk_sigma and the raw realisations .jk_real
                 (shape (nbins*nmu, njk, ncols)) to the result.

    Returns
    -------
    TwoPCFResult
    """
    binary = _find_binary(binary)
    randoms_pos, randoms_weights = _prepare_randoms(
        randoms_pos, randoms_weights, box, 'RR')
    with tempfile.TemporaryDirectory() as tmp:
        gal = os.path.join(tmp, 'data.gal')
        out = os.path.join(tmp, 'result.out')
        _write_catalog(gal, positions, weights)
        args = [binary, '-gal', gal,
                '-rmin', str(rmin), '-rmax', str(rmax),
                '-nbins', str(nbins), '-nmu', str(nmu),
                '-out', out, '-2pcf']
        if njk:
            args += ['-njk', str(njk)]
        if randoms_pos is not None:
            ran = os.path.join(tmp, 'rand.ran')
            _write_catalog(ran, randoms_pos, randoms_weights)
            args += ['-ran', ran]
        if box is not None:
            args += ['-box', _box_arg(box)]
        stdout = _run(args, tmp)
        result = TwoPCFResult(_load(out, stdout))
        # Legendre multipoles are written whenever a per-pair mu exists
        # (periodic box: plane-parallel z; survey mode: nmu > 1).
        if os.path.exists(out + '.mult'):
            mult = _load(out + '.mult', stdout)
            result.xi0 = mult[:, 2]
            result.xi2 = mult[:, 3]
            result.xi4 = mult[:, 4]
        if njk:
            real, err = _load_jk(out, stdout, njk)
            result.njk = njk
            result.jk_real = real
            result.xi_jk_mean = err[:, -2]
            result.xi_jk_sigma = err[:, -1]
        return result


def compute_3pcf(positions, weights, randoms_pos=None, randoms_weights=None,
                 rmin=1.0, rmax=30.0, nbins=6, nmu=1, box=None, binary=None,
                 njk=0):
    """Compute the 3-point correlation function (all triangle configurations).

    Parameters
    ----------
    positions : array (N, 3)
    weights   : array (N,)
    randoms_pos, randoms_weights : optional random catalog
    rmin, rmax : float
    nbins      : int   — number of radial bins (n_configs ≈ nbins³/6)
    nmu        : int   — number of mu bins (1 = isotropic)
    box        : float or (3,) — periodic box side length(s).  Separations use
                 the minimum image; if no randoms are given, RRR is computed
                 analytically and no random catalog is needed.  Requires
                 rmax <= L/4 and nmu = 1.
    binary     : str or Path
    njk        : int — number of delete-one angular jackknife regions
                 (0 = off; requires randoms, incompatible with box).  Adds
                 .zeta_jk_mean, .zeta_jk_sigma and .jk_real to the result.

    Returns
    -------
    ThreePCFResult
    """
    binary = _find_binary(binary)
    randoms_pos, randoms_weights = _prepare_randoms(
        randoms_pos, randoms_weights, box, 'RRR')
    with tempfile.TemporaryDirectory() as tmp:
        gal = os.path.join(tmp, 'data.gal')
        out = os.path.join(tmp, 'result.out')
        _write_catalog(gal, positions, weights)
        args = [binary, '-gal', gal,
                '-rmin', str(rmin), '-rmax', str(rmax),
                '-nbins', str(nbins), '-nmu', str(nmu),
                '-out', out, '-3pcf']
        if njk:
            args += ['-njk', str(njk)]
        if randoms_pos is not None:
            ran = os.path.join(tmp, 'rand.ran')
            _write_catalog(ran, randoms_pos, randoms_weights)
            args += ['-ran', ran]
        if box is not None:
            args += ['-box', _box_arg(box)]
        stdout = _run(args, tmp)
        result = ThreePCFResult(_load(out, stdout))
        if njk:
            real, err = _load_jk(out, stdout, njk)
            result.njk = njk
            result.jk_real = real
            result.zeta_jk_mean = err[:, -2]
            result.zeta_jk_sigma = err[:, -1]
        return result


def compute_4pcf(positions, weights, randoms_pos=None, randoms_weights=None,
                 rmin=1.0, rmax=30.0, nbins=3, parity=False, box=None,
                 binary=None, ntheta=None, nphi=None, exact_parity=False,
                 njk=0):
    """Compute the 4-point correlation function.

    Parameters
    ----------
    positions : array (N, 3)
    weights   : array (N,)
    randoms_pos, randoms_weights : optional random catalog
    rmin, rmax : float
    nbins      : int   — number of radial bins (n_configs grows as ~nbins⁶/48)
    parity     : bool  — if True, decompose into parity-even and parity-odd channels
    box        : float or (3,) — periodic box side length(s).  Separations use
                 the minimum image; if no randoms are given, RRRR is computed
                 analytically and no random catalog is needed.  Requires
                 rmax <= L/4.
    binary     : str or Path
    ntheta, nphi : int, optional — direction-pixel grid for the parity channel
        (default 4 x 16; the product must be <= 32767).  Finer grids reduce the
        chirality attenuation of near-degenerate tetrahedra.
    exact_parity : bool — compute the parity sign from the exact galaxy
        positions instead of pixelized directions.  No attenuation and no
        discarded tetrahedra; ignores ntheta/nphi.
    njk : int — number of delete-one angular jackknife regions (0 = off;
        requires randoms, incompatible with box).  Adds .zeta_jk_mean,
        .zeta_jk_sigma, .zeta_conn_jk_mean, .zeta_conn_jk_sigma (and with
        parity=True .zeta_odd_jk_mean, .zeta_odd_jk_sigma) plus the raw
        realisations .jk_real to the result.

    Returns
    -------
    FourPCFResult
        .zeta_conn   — connected 4PCF  (parity=False)
        .zeta_even   — parity-even 4PCF (parity=True)
        .zeta_odd    — parity-odd 4PCF  (parity=True)
    """
    binary = _find_binary(binary)
    randoms_pos, randoms_weights = _prepare_randoms(
        randoms_pos, randoms_weights, box, 'RRRR')
    with tempfile.TemporaryDirectory() as tmp:
        gal = os.path.join(tmp, 'data.gal')
        out = os.path.join(tmp, 'result.out')
        _write_catalog(gal, positions, weights)
        flag = '-4pcfp' if parity else '-4pcf'
        args = [binary, '-gal', gal,
                '-rmin', str(rmin), '-rmax', str(rmax),
                '-nbins', str(nbins),
                '-out', out, flag]
        if njk:
            args += ['-njk', str(njk)]
        if parity:
            if exact_parity:
                args.append('-exactparity')
            if ntheta is not None:
                args += ['-ntheta', str(ntheta)]
            if nphi is not None:
                args += ['-nphi', str(nphi)]
        if randoms_pos is not None:
            ran = os.path.join(tmp, 'rand.ran')
            _write_catalog(ran, randoms_pos, randoms_weights)
            args += ['-ran', ran]
        if box is not None:
            args += ['-box', _box_arg(box)]
        stdout = _run(args, tmp)
        result = FourPCFResult(_load(out, stdout), parity=parity)
        if njk:
            real, err = _load_jk(out, stdout, njk)
            result.njk = njk
            result.jk_real = real
            if parity:
                result.zeta_jk_mean       = err[:, 12]
                result.zeta_jk_sigma      = err[:, 13]
                result.zeta_odd_jk_mean   = err[:, 14]
                result.zeta_odd_jk_sigma  = err[:, 15]
                result.zeta_conn_jk_mean  = err[:, 16]
                result.zeta_conn_jk_sigma = err[:, 17]
            else:
                result.zeta_jk_mean       = err[:, 12]
                result.zeta_jk_sigma      = err[:, 13]
                result.zeta_conn_jk_mean  = err[:, 14]
                result.zeta_conn_jk_sigma = err[:, 15]
        return result
