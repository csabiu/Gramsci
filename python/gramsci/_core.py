import os
import subprocess
import tempfile
import numpy as np
from pathlib import Path

_REPO_ROOT = Path(__file__).parent.parent.parent
_DEFAULT_BINARY = _REPO_ROOT / 'bin' / 'gramsci'


def _find_binary(binary):
    if binary is not None:
        return str(binary)
    env = os.environ.get('GRAMSCI_BIN')
    if env:
        return env
    if _DEFAULT_BINARY.exists():
        return str(_DEFAULT_BINARY)
    return 'gramsci'


def _write_catalog(path, positions, weights):
    data = np.column_stack([np.asarray(positions, dtype=np.float64),
                            np.asarray(weights, dtype=np.float64)])
    np.savetxt(path, data, fmt='%.10e')


def _run(args, cwd):
    result = subprocess.run(args, capture_output=True, text=True, cwd=cwd)
    if result.returncode not in (0, 1):  # Fortran STOP returns 0 or 1
        raise RuntimeError(
            f"gramsci exited with code {result.returncode}:\n"
            f"{result.stdout[-2000:]}\n{result.stderr[-500:]}"
        )
    return result.stdout


def _load(path):
    if not os.path.exists(path):
        raise FileNotFoundError(
            f"gramsci did not produce output at {path}. "
            "Check that the binary ran successfully."
        )
    return np.loadtxt(path, skiprows=1)


# ---------------------------------------------------------------------------
# Result containers
# ---------------------------------------------------------------------------

class TwoPCFResult:
    """Result of compute_2pcf.

    Attributes
    ----------
    r : midpoint of each radial bin
    r_min, r_max : bin edges
    NN, RR : weighted pair counts
    xi : correlation function estimate
    """
    def __init__(self, data):
        # columns: r_min r_max mu_min mu_max NN RR xi
        self.r_min = data[:, 0]
        self.r_max = data[:, 1]
        self.NN    = data[:, 4]
        self.RR    = data[:, 5]
        self.xi    = data[:, 6]
        self.r     = 0.5 * (self.r_min + self.r_max)

    def __repr__(self):
        return f"TwoPCFResult(nbins={len(self.r)})"


class ThreePCFResult:
    """Result of compute_3pcf (all triangle configurations).

    Attributes
    ----------
    r1, r2, r3 : midpoints of triangle side bins
    r{1,2,3}_min, r{1,2,3}_max : bin edges
    NNN, RRR : weighted triplet counts
    zeta : 3PCF estimate  zeta = NNN / RRR
    n_configs : number of distinct (r1,r2,r3) configurations
    """
    def __init__(self, data):
        self.r1_min = data[:, 0]; self.r1_max = data[:, 1]
        self.r2_min = data[:, 2]; self.r2_max = data[:, 3]
        self.r3_min = data[:, 4]; self.r3_max = data[:, 5]
        self.NNN    = data[:, 6]
        self.RRR    = data[:, 7]
        self.zeta   = data[:, 8]
        self.r1     = 0.5 * (self.r1_min + self.r1_max)
        self.r2     = 0.5 * (self.r2_min + self.r2_max)
        self.r3     = 0.5 * (self.r3_min + self.r3_max)
        self.n_configs = len(self.zeta)

    def __repr__(self):
        return f"ThreePCFResult(n_configs={self.n_configs})"


class FourPCFResult:
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

    def __repr__(self):
        kind = 'parity' if self.parity else 'connected'
        return f"FourPCFResult(n_configs={self.n_configs}, kind={kind!r})"


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------

def compute_2pcf(positions, weights, randoms_pos=None, randoms_weights=None,
                 rmin=1.0, rmax=30.0, nbins=10, nmu=1, binary=None):
    """Compute the 2-point correlation function.

    Parameters
    ----------
    positions : array (N, 3)  — x, y, z coordinates
    weights   : array (N,)    — per-particle weights
    randoms_pos, randoms_weights : optional random catalog for RR estimation
    rmin, rmax : float        — radial range
    nbins      : int          — number of radial bins
    nmu        : int          — number of mu bins (1 = isotropic)
    binary     : str or Path  — path to gramsci binary (auto-detected if None)

    Returns
    -------
    TwoPCFResult
    """
    binary = _find_binary(binary)
    with tempfile.TemporaryDirectory() as tmp:
        gal = os.path.join(tmp, 'data.gal')
        out = os.path.join(tmp, 'result.out')
        _write_catalog(gal, positions, weights)
        args = [binary, '-gal', gal,
                '-rmin', str(rmin), '-rmax', str(rmax),
                '-nbins', str(nbins), '-nmu', str(nmu),
                '-out', out, '-2pcf']
        if randoms_pos is not None:
            ran = os.path.join(tmp, 'rand.ran')
            _write_catalog(ran, randoms_pos, randoms_weights)
            args += ['-ran', ran]
        _run(args, tmp)
        return TwoPCFResult(_load(out))


def compute_3pcf(positions, weights, randoms_pos=None, randoms_weights=None,
                 rmin=1.0, rmax=30.0, nbins=6, nmu=1, binary=None):
    """Compute the 3-point correlation function (all triangle configurations).

    Parameters
    ----------
    positions : array (N, 3)
    weights   : array (N,)
    randoms_pos, randoms_weights : optional random catalog
    rmin, rmax : float
    nbins      : int   — number of radial bins (n_configs ≈ nbins³/6)
    nmu        : int   — number of mu bins (1 = isotropic)
    binary     : str or Path

    Returns
    -------
    ThreePCFResult
    """
    binary = _find_binary(binary)
    with tempfile.TemporaryDirectory() as tmp:
        gal = os.path.join(tmp, 'data.gal')
        out = os.path.join(tmp, 'result.out')
        _write_catalog(gal, positions, weights)
        args = [binary, '-gal', gal,
                '-rmin', str(rmin), '-rmax', str(rmax),
                '-nbins', str(nbins), '-nmu', str(nmu),
                '-out', out, '-3pcf']
        if randoms_pos is not None:
            ran = os.path.join(tmp, 'rand.ran')
            _write_catalog(ran, randoms_pos, randoms_weights)
            args += ['-ran', ran]
        _run(args, tmp)
        return ThreePCFResult(_load(out))


def compute_4pcf(positions, weights, randoms_pos=None, randoms_weights=None,
                 rmin=1.0, rmax=30.0, nbins=3, parity=False, binary=None):
    """Compute the 4-point correlation function.

    Parameters
    ----------
    positions : array (N, 3)
    weights   : array (N,)
    randoms_pos, randoms_weights : optional random catalog
    rmin, rmax : float
    nbins      : int   — number of radial bins (n_configs grows as ~nbins⁶/48)
    parity     : bool  — if True, decompose into parity-even and parity-odd channels
    binary     : str or Path

    Returns
    -------
    FourPCFResult
        .zeta_conn   — connected 4PCF  (parity=False)
        .zeta_even   — parity-even 4PCF (parity=True)
        .zeta_odd    — parity-odd 4PCF  (parity=True)
    """
    binary = _find_binary(binary)
    with tempfile.TemporaryDirectory() as tmp:
        gal = os.path.join(tmp, 'data.gal')
        out = os.path.join(tmp, 'result.out')
        _write_catalog(gal, positions, weights)
        flag = '-4pcfp' if parity else '-4pcf'
        args = [binary, '-gal', gal,
                '-rmin', str(rmin), '-rmax', str(rmax),
                '-nbins', str(nbins),
                '-out', out, flag]
        if randoms_pos is not None:
            ran = os.path.join(tmp, 'rand.ran')
            _write_catalog(ran, randoms_pos, randoms_weights)
            args += ['-ran', ran]
        _run(args, tmp)
        return FourPCFResult(_load(out), parity=parity)
