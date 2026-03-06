"""countdata._core — ctypes interface to the C beta-binomial library."""

import ctypes
import multiprocessing
import os
import platform

import numpy as np

# ── Load shared library ────────────────────────────────────────────────────

def _find_lib():
    pkg = os.path.dirname(os.path.abspath(__file__))
    names = (["_countdata.dylib", "_countdata.so", "libcountdata.so"]
             if platform.system() == "Darwin"
             else ["_countdata.so", "libcountdata.so"])
    for n in names:
        for d in (pkg, os.path.dirname(pkg), os.getcwd()):
            p = os.path.join(d, n)
            if os.path.isfile(p):
                return p
    raise FileNotFoundError(
        f"_countdata.so not found in {pkg}. "
        "Reinstall: pip install --force-reinstall pycountdata"
    )

_lib = ctypes.CDLL(_find_lib())

_P_INT = ctypes.POINTER(ctypes.c_int)
_P_DBL = ctypes.POINTER(ctypes.c_double)

_lib.bb.restype = None
_lib.bb.argtypes = [_P_INT, _P_DBL, _P_DBL, _P_INT, _P_INT, _P_INT, _P_DBL, _P_INT, _P_DBL]

_lib.ibb.restype = None
_lib.ibb.argtypes = [_P_INT, _P_DBL, _P_DBL, _P_DBL, _P_DBL, _P_INT, _P_DBL, _P_INT, _P_DBL, _P_DBL]

_lib.bbCores.restype = None
_lib.bbCores.argtypes = [_P_INT]

# ── Thread helpers ──────────────────────────────────────────────────────────

def _ncores():
    n = ctypes.c_int(0)
    _lib.bbCores(ctypes.byref(n))
    return n.value or multiprocessing.cpu_count()

def _threads(n, K):
    if n == 1:
        return 1
    nc = _ncores()
    t = min(n, nc) if n > 0 else max(1, nc + n)
    return min(t, K)

def _ptr(arr):
    return arr.ctypes.data_as(_P_DBL)

def _iptr(arr):
    return arr.ctypes.data_as(_P_INT)

# ── bb_test ─────────────────────────────────────────────────────────────────

def bb_test(x, tx, group, alternative="two.sided", n_threads=-1, verbose=True):
    """Beta-binomial test for differential count data (unpaired).

    Parameters
    ----------
    x : array_like, shape (K, N) or (N,)
        Count matrix. K features (genes), N samples.
    tx : array_like, shape (K, N) or (N,)
        Total counts per sample. Must be positive.
    group : array_like, shape (N,)
        Group labels. At least 2 groups.
    alternative : {'two.sided', 'less', 'greater'}
        One-sided only for exactly 2 groups.
    n_threads : int
        Number of threads. -1 = all cores.
    verbose : bool
        Print progress to stdout.

    Returns
    -------
    dict  —  ``{'p.value': ndarray}`` of length K.
    """
    x, tx, group = _validate(x, tx, group)
    K, N = x.shape

    gvals = list(dict.fromkeys(group))
    M = len(gvals)
    if M < 2:
        raise ValueError("Need at least 2 groups.")
    if alternative not in ("two.sided", "less", "greater"):
        raise ValueError("alternative must be 'two.sided', 'less', or 'greater'")
    if alternative != "two.sided" and M != 2:
        raise ValueError("One-sided test for two groups only.")

    g_idx = [np.where(group == g)[0] for g in gvals]
    g_size = np.array([len(g) for g in g_idx], dtype=np.int32)
    g_ind = np.zeros(M, dtype=np.int32)
    for i in range(1, M):
        g_ind[i] = g_ind[i - 1] + g_size[i - 1]

    # flatten in group order (matches R's .C call layout)
    a = np.empty(K * N, dtype=np.float64)
    ta = np.empty(K * N, dtype=np.float64)
    p = 0
    for k in range(K):
        for m in range(M):
            s = g_size[m]
            a[p:p + s] = x[k, g_idx[m]]
            ta[p:p + s] = tx[k, g_idx[m]]
            p += s

    nt = _threads(n_threads, K)
    mem = np.zeros(nt * (2 * N + M), dtype=np.float64)
    mem[0] = 1.0  # theta_equal
    mem[1] = {"two.sided": 0.0, "less": -1.0, "greater": 1.0}[alternative]

    pval = np.zeros(K, dtype=np.float64)
    _lib.bb(ctypes.c_int(K), _ptr(a), _ptr(ta), ctypes.c_int(M),
            _iptr(g_size), _iptr(g_ind), _ptr(mem),
            ctypes.c_int(nt if verbose else -nt), _ptr(pval))
    return {"p.value": pval}


# ── ibb_test ────────────────────────────────────────────────────────────────

def ibb_test(x, tx, group, alternative="two.sided", n_threads=-1,
             BIG=1e4, verbose=True):
    """Inverted beta-binomial test for paired count data.

    Parameters
    ----------
    x : array_like, shape (K, N) or (N,)
        Count matrix.
    tx : array_like, shape (K, N) or (N,)
        Total counts. Must be positive.
    group : array_like, shape (N,)
        Exactly 2 groups with equal sample sizes.
    alternative : {'two.sided', 'less', 'greater'}
    n_threads : int
        -1 = all cores.
    BIG : float
        Fold-change cap for zero-count edge cases.
    verbose : bool

    Returns
    -------
    dict  —  ``{'p.value': ndarray, 'fc': ndarray}``.
    """
    x, tx, group = _validate(x, tx, group)
    K, N = x.shape

    gvals = list(dict.fromkeys(group))
    if len(gvals) != 2:
        raise ValueError("Paired test requires exactly 2 groups.")
    if alternative not in ("two.sided", "less", "greater"):
        raise ValueError("alternative must be 'two.sided', 'less', or 'greater'")

    ia = np.where(group == gvals[0])[0]
    ib = np.where(group == gvals[1])[0]
    if len(ia) != len(ib):
        raise ValueError("Groups must have equal sizes for paired test.")

    Np = len(ia)
    a = np.empty(K * Np, dtype=np.float64)
    b = np.empty(K * Np, dtype=np.float64)
    ta = np.empty(K * Np, dtype=np.float64)
    tb = np.empty(K * Np, dtype=np.float64)
    for i in range(K):
        s, e = i * Np, (i + 1) * Np
        a[s:e], b[s:e] = x[i, ia], x[i, ib]
        ta[s:e], tb[s:e] = tx[i, ia], tx[i, ib]

    Q = 2 ** 15
    nt = _threads(n_threads, K)
    mem = np.zeros(nt * (4 * Np + Np * Q + 5 * Q) + 3 * Q, dtype=np.float64)
    mem[0] = 5.0  # lower_bound
    mem[1] = {"two.sided": 0.0, "less": -1.0, "greater": 1.0}[alternative]

    pval = np.zeros(K, dtype=np.float64)
    fc = np.zeros(K, dtype=np.float64)

    _lib.ibb(ctypes.c_int(K), _ptr(a), _ptr(b), _ptr(ta), _ptr(tb),
             ctypes.c_int(Np), _ptr(mem),
             ctypes.c_int(nt if verbose else -nt), _ptr(pval), _ptr(fc))

    # post-process fold change (matches R)
    fc[fc < 1.0] = -1.0 / fc[fc < 1.0]
    if K > 1:
        t1 = x[:, ia].sum(axis=1) if Np > 1 else x[:, ia].ravel()
        t2 = x[:, ib].sum(axis=1) if Np > 1 else x[:, ib].ravel()
        fc[t1 == 0] = fc[t1 == 0] * 0 + BIG
        fc[t2 == 0] = fc[t2 == 0] * 0 - BIG
        fc[(t1 == 0) & (t2 == 0)] = 1.0

    return {"p.value": pval, "fc": fc}


# ── Utilities ───────────────────────────────────────────────────────────────

def normalize(d):
    """Normalize count matrix so all columns have the same total."""
    d = np.asarray(d, dtype=np.float64)
    s = d.sum(axis=0)
    return d / (s / s.mean())[np.newaxis, :]


def fold_change(d1, d2, BIG=1e4):
    """Per-row fold change. Negative = d1 > d2, positive = d2 > d1."""
    d1 = np.atleast_2d(np.asarray(d1, dtype=np.float64))
    d2 = np.atleast_2d(np.asarray(d2, dtype=np.float64))
    v1, v2 = d1.mean(axis=1), d2.mean(axis=1)
    fc = np.ones(len(v1))
    for i in range(len(v1)):
        if v1[i] == 0:
            fc[i] = BIG if v2[i] else 1.0
        elif v2[i] == 0:
            fc[i] = -BIG
        elif v1[i] > v2[i]:
            fc[i] = -v1[i] / v2[i]
        elif v1[i] < v2[i]:
            fc[i] = v2[i] / v1[i]
    return fc


# ── Input validation ────────────────────────────────────────────────────────

def _validate(x, tx, group):
    x = np.asarray(x, dtype=np.float64)
    tx = np.asarray(tx, dtype=np.float64)
    group = np.asarray(group)
    if np.any(tx <= 0):
        raise ValueError("tx must be positive.")
    if np.any(x < 0):
        raise ValueError("x must be non-negative.")
    if x.ndim == 1:
        x = x.reshape(1, -1)
    K, N = x.shape
    if tx.ndim == 1:
        tx = np.tile(tx, (K, 1))
    if tx.shape[0] == 1 and K > 1:
        tx = np.tile(tx, (K, 1))
    if len(group) != N:
        raise ValueError("length of 'group' must equal number of columns of 'x'")
    if tx.shape[1] != N:
        raise ValueError("columns of 'tx' must equal columns of 'x'")
    return x, tx, group