"""Fast critical-plane search for multiaxial fatigue on FE stress fields (solver-agnostic: NGSolve, CalculiX, ...).

    from cplane import Findley, search_prop, search_path
    res = search_prop(Sa, Sm, Findley(k=0.3))        # constant amplitude, proportional (one load channel)
    res = search_path(H, Findley(k=0.3))             # one cycle of a general (non-proportional) path

Stress layout: Voigt TENSOR components [xx, yy, zz, xy, yz, zx] (no factor 2; CalculiX .frd order).
Abaqus / ccx .dat order is [11, 22, 33, 12, 13, 23] -> use `from_abaqus`.

Theory.  On a plane with unit normal n and in-plane basis (u, v): sigma_n = n.S.n, tau = (u.S.n, v.S.n).
    proportional constant amplitude, S(t) = Sm + Sa f(t), f spanning [-1, 1]:  tau_a = |tau(Sa)|, sn_a = |sn(Sa)|, sn_m = sn(Sm)
    general path S(t_i):  tau_a = radius of the minimum circumscribed circle (MCC) of the shear path,
                          sn_a = (max - min)/2, sn_m = (max + min)/2.
The criterion maps (tau_a, sn_a, sn_m, |n.d|) to a scalar damage objective; the critical plane maximises it, so
build-direction-dependent strength (|n.d|, d = build direction) changes the plane choice and not only the value.

Numerics (why it is fast).
  1. Everything per point runs in one Numba kernel, parallel over points, with no temporaries: traction t = S n,
     sigma_n = n.t, tau = t - sigma_n n (proportional: 2 x 9 multiply-adds per plane).
  2. Coarse Fibonacci hemisphere (default 800 normals prop / 400 path, ~5-7 deg) -> its local maxima on an
     8-nearest-neighbour graph -> the q best of them -> Newton on the tangent plane with a 3 x 3 finite-difference stencil (compass move
     as fallback, step halving) down to `tol` rad.
     (Findley often has several peaks within 0.1-1 % of each other, tens of degrees apart, separated by saddles
     shallower than the coarse sampling error at 400 normals; 800 resolves them, see CPLANE_SPEED.md.)
     The stencil is rotated by the golden angle every iteration (dense directions, MADS-like), so the search does not
     stall on the kinks of the MCC radius (a min-max function of n) in the non-proportional case.
  3. Hot-spot screening: a cheap per-point upper bound (principal values only) orders the points; exact searches stop
     once the bound drops below `screen` x the worst exact damage found.  Unsearched points are returned as NaN.
  4. Shear amplitude = minimum circumscribed circle (MCC; one of several definitions - longest chord, MRH, ...;
     MCC under-predicts for cross / box paths). Incremental Welzl, deterministic order (O(T) expected only for
     random order; worst case slower), points translated to the first sample for precision.
  The vectorised NumPy projections (`proj_rows`, `_quant_prop`, `_quant_path`) are kept as an independent reference
  for brute-force checks.
Anisotropy: strength(c) = sum_i poly[i] c^i with c = |n.d|; the objective is divided by it (must stay > 0).
"""
import numpy as np

try:
    import numba
except ImportError:                                           # pragma: no cover
    numba = None

VOIGT = ("xx", "yy", "zz", "xy", "yz", "zx")


def from_abaqus(S):
    """[11, 22, 33, 12, 13, 23] -> [xx, yy, zz, xy, yz, zx]."""
    S = np.asarray(S)
    return S[..., [0, 1, 2, 3, 5, 4]]


def to_matrix(S):
    S = np.asarray(S)
    xx, yy, zz, xy, yz, zx = np.moveaxis(S, -1, 0)
    return np.stack([np.stack([xx, xy, zx], -1), np.stack([xy, yy, yz], -1), np.stack([zx, yz, zz], -1)], -2)


def to_voigt(A):
    return np.stack([A[..., 0, 0], A[..., 1, 1], A[..., 2, 2], A[..., 0, 1], A[..., 1, 2], A[..., 2, 0]], -1)


# ------------------------------------------------------------------------------------------------- plane geometry
def fib_hemisphere(m):
    """m quasi-uniform unit normals on the upper hemisphere (n and -n are the same plane). Spacing ~ sqrt(2 pi / m)."""
    i = np.arange(m) + 0.5
    z = 1.0 - i / m
    r = np.sqrt(1.0 - z * z)
    phi = i * np.pi * (3.0 - np.sqrt(5.0))
    return np.stack([r * np.cos(phi), r * np.sin(phi), z], -1)


_GRIDS = {}


def coarse_grid(m, k=8):
    """Fibonacci hemisphere + its k-nearest-neighbour graph (|n.m|, so n and -n coincide across the equator)."""
    if m not in _GRIDS:
        N = fib_hemisphere(m)
        C = np.abs(N @ N.T)
        np.fill_diagonal(C, -1.0)
        _GRIDS[m] = (N, np.ascontiguousarray(np.argsort(-C, axis=1)[:, :k]))
    return _GRIDS[m]


def tangent_basis(n):
    a = np.zeros_like(n)
    small = np.abs(n[..., 0]) < 0.9
    a[..., 0] = small
    a[..., 1] = ~small
    u = np.cross(n, a)
    u /= np.linalg.norm(u, axis=-1, keepdims=True)
    return u, np.cross(n, u)


def _row(a, b):
    """Voigt row r with r . S = a . S . b for symmetric S."""
    return np.stack([a[..., 0] * b[..., 0], a[..., 1] * b[..., 1], a[..., 2] * b[..., 2],
                     a[..., 0] * b[..., 1] + a[..., 1] * b[..., 0], a[..., 1] * b[..., 2] + a[..., 2] * b[..., 1],
                     a[..., 2] * b[..., 0] + a[..., 0] * b[..., 2]], -1)


def proj_rows(n):
    """(..., 3, 6): rows giving (sigma_n, tau_u, tau_v) on the plane with normal n."""
    u, v = tangent_basis(n)
    return np.stack([_row(n, n), _row(u, n), _row(v, n)], -2)


# --------------------------------------------------------------------------------------------------- criteria
class Criterion:
    """Objective on plane quantities divided by an orientation-dependent strength polynomial in c = |n.d|.
    Note: Findley values can be negative (compressive mean); dividing a negative value by a larger strength makes it
    less negative, i.e. 'more damaging' - only positive values are physically meaningful with anisotropy."""
    name, code, prm = "base", -1, 0.0

    def __init__(self, build_dir=None, strength_poly=None):
        self.d = np.zeros(3) if build_dir is None else np.asarray(build_dir, float) / np.linalg.norm(build_dir)
        self.poly = np.zeros(0) if strength_poly is None else np.asarray(strength_poly, float)
        if self.poly.size and self.min_strength() <= 0:
            raise ValueError("strength polynomial must be > 0 on [0, 1]")

    def min_strength(self):
        """Exact minimum of the strength polynomial on [0, 1] (endpoints + real roots of p'), shaved by 1e-12 so the
        screening bound stays a true upper bound."""
        if not self.poly.size:
            return 1.0
        P = np.polynomial.Polynomial(self.poly)
        c = [0.0, 1.0] + [float(r.real) for r in P.deriv().roots() if abs(r.imag) < 1e-12 and 0.0 <= r.real <= 1.0]
        return float(min(P(np.array(c)))) * (1 - 1e-12)

    def strength(self, c):
        return np.polynomial.polynomial.polyval(c, self.poly) if self.poly.size else np.ones_like(c)

    def objective(self, tau_a, sn_a, sn_m, n):
        return self.raw(tau_a, sn_a, sn_m) / self.strength(np.abs(n @ self.d))

    def bound(self, TB, AB, NB):                     # upper bound from tau_a <= TB, sn_a <= AB, sn_max <= NB
        return self.raw_bound(TB, AB, NB) / self.min_strength()


class Findley(Criterion):
    """f = tau_a + k sigma_n,max   (Findley 1959; shear cracks, high-cycle).  k >= 0."""
    name, code = "findley", 0

    def __init__(self, k, **kw):
        super().__init__(**kw)
        if k < 0:
            raise ValueError("Findley k must be >= 0")
        self.prm = self.k = float(k)

    def raw(self, tau_a, sn_a, sn_m):
        return tau_a + self.k * (sn_m + sn_a)

    def raw_bound(self, TB, AB, NB):
        return TB + self.k * np.maximum(NB, 0.0)


class MWCM(Criterion):
    """Modified Woehler Curve Method (Susmel & Lazzarin 2002): plane of maximum tau_a, ties broken by the largest
    sigma_n,max (eps-weighted: the plane tilts by O(eps), so rho is biased O(eps) and tau_a O(eps^2); eps = 1e-6).
    Life via `mwcm_life`."""
    name, code, prm = "mwcm", 1, 1e-6

    def raw(self, tau_a, sn_a, sn_m):
        return tau_a + self.prm * (sn_m + sn_a)

    def raw_bound(self, TB, AB, NB):
        return TB + self.prm * np.maximum(NB, 0.0)


class SWTStress(Criterion):
    """Tensile-crack plane: sqrt(sigma_n,max * sigma_n,a) (elastic Smith-Watson-Topper, stress units).
    Suited to inter-layer (mode I) failure of printed polymers."""
    name, code = "swt", 2

    def raw(self, tau_a, sn_a, sn_m):
        return np.sqrt(np.maximum(sn_m + sn_a, 0.0) * sn_a)

    def raw_bound(self, TB, AB, NB):
        return np.sqrt(np.maximum(NB, 0.0) * AB)


class WalkerStress(Criterion):
    """Walker-type plane parameter sigma_n,max^(1 - g) * sigma_n,a^g (stress units, homogeneous of degree 1).
    g = 0: max normal stress, g = 0.5: SWT-stress, g = 1: normal-stress amplitude (mean-stress insensitive).
    Zero when sigma_n,max <= 0 (closed crack plane)."""
    name, code = "walker", 3

    def __init__(self, gamma, **kw):
        super().__init__(**kw)
        if not 0.0 <= gamma <= 1.0:
            raise ValueError("Walker gamma must be in [0, 1]")
        self.prm = self.gamma = float(gamma)

    def raw(self, tau_a, sn_a, sn_m):
        smax = sn_m + sn_a
        with np.errstate(divide="ignore", invalid="ignore"):
            v = np.maximum(smax, 0.0) ** (1 - self.gamma) * sn_a ** self.gamma
        return np.where(smax > 0, v, 0.0)

    def raw_bound(self, TB, AB, NB):
        return np.maximum(NB, 0.0) ** (1 - self.gamma) * AB ** self.gamma


class FatemiSocie(Criterion):
    """Fatemi-Socie in stress units for linear elasticity: tau_a (1 + k sigma_n,max / sigma_y).
    Engineering shear-strain amplitude = tau_a / G on every plane, so FS = value / G.  Here the plane maximises the
    whole parameter (the 'extended' FS of Chiocca et al.); the classic choice (plane of max shear amplitude) is
    MWCM's plane - evaluate `raw` on its tau_a, sn_a, sn_m."""
    name, code = "fs", 4

    def __init__(self, k, sy, **kw):
        super().__init__(**kw)
        if k < 0 or sy <= 0:
            raise ValueError("FS needs k >= 0, sigma_y > 0")
        self.k, self.sy = float(k), float(sy)
        self.prm = self.k / self.sy

    def raw(self, tau_a, sn_a, sn_m):
        return tau_a * (1.0 + self.prm * (sn_m + sn_a))

    def raw_bound(self, TB, AB, NB):
        return TB * (1.0 + self.prm * np.maximum(NB, 0.0))


def mwcm_life(tau_a, sn_a, sn_m, tau_A_of_rho, k_of_rho, NA=2e6, m=1.0, rho_lim=None):
    """N = NA (tau_A,ref(rho) / tau_a)^k(rho),  rho_eff = (m sn_m + sn_a) / tau_a (capped at rho_lim)."""
    rho = (m * sn_m + sn_a) / np.maximum(tau_a, 1e-300)
    if rho_lim is not None:
        rho = np.minimum(rho, rho_lim)
    return NA * (tau_A_of_rho(rho) / tau_a) ** k_of_rho(rho), rho


# ---------------------------------------------------------------------------------------- minimum circumscribed circle
if numba is not None:
    @numba.njit(cache=True, fastmath=False)
    def _mcc(x, y):
        T = x.size
        cx, cy, r2 = x[0], y[0], 0.0
        for i in range(1, T):
            if (x[i] - cx) ** 2 + (y[i] - cy) ** 2 <= r2 * (1 + 1e-12) + 1e-300:
                continue
            cx, cy, r2 = x[i], y[i], 0.0
            for j in range(i):
                if (x[j] - cx) ** 2 + (y[j] - cy) ** 2 <= r2 * (1 + 1e-12) + 1e-300:
                    continue
                cx, cy = 0.5 * (x[i] + x[j]), 0.5 * (y[i] + y[j])
                r2 = (x[i] - cx) ** 2 + (y[i] - cy) ** 2
                for k in range(j):
                    if (x[k] - cx) ** 2 + (y[k] - cy) ** 2 <= r2 * (1 + 1e-12) + 1e-300:
                        continue
                    ax, ay = x[i], y[i]
                    bx, by = x[j] - ax, y[j] - ay
                    qx, qy = x[k] - ax, y[k] - ay
                    d = 2.0 * (bx * qy - by * qx)
                    if abs(d) < 1e-14 * (bx * bx + by * by + qx * qx + qy * qy) + 1e-300:   # collinear: farthest pair
                        best, ux, uy, vx, vy = -1.0, 0.0, 0.0, 0.0, 0.0
                        for (px, py, sx, sy) in ((x[i], y[i], x[j], y[j]), (x[i], y[i], x[k], y[k]),
                                                 (x[j], y[j], x[k], y[k])):
                            dd = (px - sx) ** 2 + (py - sy) ** 2
                            if dd > best:
                                best, ux, uy, vx, vy = dd, px, py, sx, sy
                        cx, cy = 0.5 * (ux + vx), 0.5 * (uy + vy)
                        r2 = 0.25 * best
                    else:
                        b2, q2 = bx * bx + by * by, qx * qx + qy * qy
                        ox = (qy * b2 - by * q2) / d
                        oy = (bx * q2 - qx * b2) / d
                        cx, cy = ax + ox, ay + oy
                        r2 = ox * ox + oy * oy
        return np.sqrt(r2)

    @numba.njit(parallel=True, cache=True)
    def _mcc_batch(U, V):
        R = U.shape[0]
        out = np.empty(R)
        for r in numba.prange(R):
            out[r] = _mcc(U[r], V[r])
        return out


def mcc_radius(U, V):
    """Minimum circumscribed circle radius of the 2D paths (U[..., t], V[..., t])."""
    if numba is None:                                          # pragma: no cover
        raise ImportError("numba is required for non-proportional paths")
    sh = U.shape[:-1]
    U = np.asarray(U, float).reshape(-1, U.shape[-1])
    V = np.asarray(V, float).reshape(-1, V.shape[-1])
    U = np.ascontiguousarray(U - U[:, :1])                    # translate: precision at large offsets
    V = np.ascontiguousarray(V - V[:, :1])
    return _mcc_batch(U, V).reshape(sh)


# ------------------------------------------------------------------------------------------------ plane quantities
def _quant_prop(W, Sa, Sm):
    """W (..., 3, 6) per point or shared (M, 3, 6). Sa, Sm (P, 6). Returns tau_a, sn_a, sn_m with shape (P, ...)."""
    if W.ndim == 3:                                            # shared normals: one GEMM
        M = W.shape[0]
        X = (Sa @ W.reshape(-1, 6).T).reshape(len(Sa), M, 3)
        snm = Sm @ W[:, 0, :].T
    else:                                                      # per-point normals (P, ..., 3, 6)
        X = np.einsum("p...cj,pj->p...c", W, Sa)
        snm = np.einsum("p...j,pj->p...", W[..., 0, :], Sm)
    return np.hypot(X[..., 1], X[..., 2]), np.abs(X[..., 0]), snm


def _quant_path(W, H):
    """H (P, T, 6) one cycle. Returns tau_a (MCC), sn_a, sn_m."""
    if W.ndim == 3:
        M = W.shape[0]
        Y = (H.reshape(-1, 6) @ W.reshape(-1, 6).T).reshape(H.shape[0], H.shape[1], M, 3)
        Y = np.moveaxis(Y, 1, -2)                              # (P, M, T, 3)
    else:
        Y = np.einsum("p...cj,ptj->p...tc", W, H)
    sn = Y[..., 0]
    hi, lo = sn.max(-1), sn.min(-1)
    return mcc_radius(Y[..., 1], Y[..., 2]), 0.5 * (hi - lo), 0.5 * (hi + lo)


# ----------------------------------------------------------------------------------------------- Numba search core
if numba is not None:
    @numba.njit(cache=True, inline="always")
    def _obj(code, prm, ta, sa, sm):
        if code == 0:
            return ta + prm * (sm + sa)
        if code == 1:
            return ta + prm * (sm + sa)
        if code == 2:
            x = (sm + sa) * sa
            return np.sqrt(x) if sm + sa > 0.0 else 0.0
        if code == 4:                                                   # Fatemi-Socie, prm = k / sigma_y
            return ta * (1.0 + prm * (sm + sa))
        smax = sm + sa                                                  # code 3: Walker, prm = gamma
        return smax ** (1.0 - prm) * sa ** prm if smax > 0.0 else 0.0

    @numba.njit(cache=True, inline="always")
    def _strength_at(poly, d, n0, n1, n2):
        if poly.size == 0:
            return 1.0
        c = abs(n0 * d[0] + n1 * d[1] + n2 * d[2])
        s, ck = 0.0, 1.0
        for a in poly:
            s += a * ck
            ck *= c
        return s

    @numba.njit(cache=True, inline="always")
    def _basis(n0, n1, n2):
        if abs(n0) < 0.9:                                      # u = n x e_x
            ux, uy, uz = 0.0, n2, -n1
        else:                                                  # u = n x e_y
            ux, uy, uz = -n2, 0.0, n0
        r = np.sqrt(ux * ux + uy * uy + uz * uz)
        ux, uy, uz = ux / r, uy / r, uz / r
        return ux, uy, uz, n1 * uz - n2 * uy, n2 * ux - n0 * uz, n0 * uy - n1 * ux

    @numba.njit(cache=True, inline="always")
    def _sn_tau(S, n0, n1, n2):
        tx = S[0] * n0 + S[3] * n1 + S[5] * n2
        ty = S[3] * n0 + S[1] * n1 + S[4] * n2
        tz = S[5] * n0 + S[4] * n1 + S[2] * n2
        sn = tx * n0 + ty * n1 + tz * n2
        return sn, tx - sn * n0, ty - sn * n1, tz - sn * n2

    @numba.njit(cache=True)
    def _eval(mode, A, B, n0, n1, n2, code, prm, d, poly, bu, bv, q3):
        """mode 0: A[0] = amplitude, B[0] = mean;  mode 1: A = path (T, 6).  q3 <- (tau_a, sn_a, sn_m)."""
        if mode == 0:
            sn, rx, ry, rz = _sn_tau(A[0], n0, n1, n2)
            ta, sa = np.sqrt(rx * rx + ry * ry + rz * rz), abs(sn)
            sm = _sn_tau(B[0], n0, n1, n2)[0]
        else:
            ux, uy, uz, vx, vy, vz = _basis(n0, n1, n2)
            hi, lo = -np.inf, np.inf
            for t in range(A.shape[0]):
                sn, rx, ry, rz = _sn_tau(A[t], n0, n1, n2)
                bu[t] = rx * ux + ry * uy + rz * uz
                bv[t] = rx * vx + ry * vy + rz * vz
                hi, lo = max(hi, sn), min(lo, sn)
            u0, v0 = bu[0], bv[0]
            for t in range(A.shape[0]):                       # translate: MCC precision at large offsets
                bu[t] -= u0
                bv[t] -= v0
            ta, sa, sm = _mcc(bu, bv), 0.5 * (hi - lo), 0.5 * (hi + lo)
        q3[0], q3[1], q3[2] = ta, sa, sm
        return _obj(code, prm, ta, sa, sm) / _strength_at(poly, d, n0, n1, n2)

    @numba.njit(cache=True)
    def _refine(mode, Ap, Bp, f, n0, n1, n2, h, tol, code, prm, d, poly, bu, bv, q3):
        """Ascent on the tangent plane: 3 x 3 stencil (rotated by the golden angle every iteration: dense directions,
        no stall on kinks), finite-difference Newton step when the stencil Hessian is negative definite, compass move
        otherwise, step halving on failure."""
        fs = np.empty(9)
        it = 0
        while h > tol and it < 200:
            it += 1
            ux, uy, uz, vx, vy, vz = _basis(n0, n1, n2)
            cr, sr = np.cos(2.399963 * it), np.sin(2.399963 * it)
            ux, uy, uz, vx, vy, vz = (cr * ux + sr * vx, cr * uy + sr * vy, cr * uz + sr * vz,
                                      cr * vx - sr * ux, cr * vy - sr * uy, cr * vz - sr * uz)
            bf, b0, b1, b2 = f, n0, n1, n2
            k = 0
            for a in (-1.0, 0.0, 1.0):
                for b in (-1.0, 0.0, 1.0):
                    if a == 0.0 and b == 0.0:
                        fs[k] = f
                        k += 1
                        continue
                    c0 = n0 + h * (a * ux + b * vx)
                    c1 = n1 + h * (a * uy + b * vy)
                    c2 = n2 + h * (a * uz + b * vz)
                    r = np.sqrt(c0 * c0 + c1 * c1 + c2 * c2)
                    c0, c1, c2 = c0 / r, c1 / r, c2 / r
                    fc = _eval(mode, Ap, Bp, c0, c1, c2, code, prm, d, poly, bu, bv, q3)
                    fs[k] = fc
                    k += 1
                    if fc > bf:
                        bf, b0, b1, b2 = fc, c0, c1, c2
            ga = (fs[7] - fs[1]) / (2 * h)                     # fs index = 3 (a + 1) + (b + 1)
            gb = (fs[5] - fs[3]) / (2 * h)
            haa = (fs[7] - 2 * f + fs[1]) / (h * h)
            hbb = (fs[5] - 2 * f + fs[3]) / (h * h)
            hab = (fs[8] - fs[6] - fs[2] + fs[0]) / (4 * h * h)
            det = haa * hbb - hab * hab
            newton = False
            step = 0.0
            if haa < 0.0 and det > 0.0:
                sa_ = -(hbb * ga - hab * gb) / det
                sb_ = -(haa * gb - hab * ga) / det
                step = np.sqrt(sa_ * sa_ + sb_ * sb_)
                if step > 2.0 * h:
                    sa_, sb_, step = sa_ * 2.0 * h / step, sb_ * 2.0 * h / step, 2.0 * h
                c0 = n0 + sa_ * ux + sb_ * vx
                c1 = n1 + sa_ * uy + sb_ * vy
                c2 = n2 + sa_ * uz + sb_ * vz
                r = np.sqrt(c0 * c0 + c1 * c1 + c2 * c2)
                c0, c1, c2 = c0 / r, c1 / r, c2 / r
                fc = _eval(mode, Ap, Bp, c0, c1, c2, code, prm, d, poly, bu, bv, q3)
                if fc >= bf:
                    bf, b0, b1, b2 = fc, c0, c1, c2
                    newton = True
            if bf > f:
                f, n0, n1, n2 = bf, b0, b1, b2
                if newton:
                    h = min(h, max(step, 0.25 * tol))
            else:
                h *= 0.5
        return f, n0, n1, n2

    @numba.njit(cache=True)
    def _point(mode, Ap, Bp, N, NB, q, h0, tol, code, prm, d, poly, out_v, out_q, out_n, p, T, M):
        bu, bv, q3 = np.empty(T), np.empty(T), np.empty(3)
        vals = np.empty(M)
        locs = np.empty(M, dtype=np.int64)
        for m in range(M):
            vals[m] = _eval(mode, Ap, Bp, N[m, 0], N[m, 1], N[m, 2], code, prm, d, poly, bu, bv, q3)
        nloc = 0                                               # coarse local maxima (antipodal-aware graph)
        for m in range(M):
            ok = True
            for k in range(NB.shape[1]):
                if vals[NB[m, k]] > vals[m]:
                    ok = False
                    break
            if ok:
                locs[nloc] = m
                nloc += 1
        lv = np.empty(nloc)
        for j in range(nloc):
            lv[j] = -vals[locs[j]]
        ordr = np.argsort(lv)
        f, n0, n1, n2 = -np.inf, 0.0, 0.0, 1.0
        for j in range(min(q, nloc)):
            i = locs[ordr[j]]
            fj, m0, m1, m2 = _refine(mode, Ap, Bp, vals[i], N[i, 0], N[i, 1], N[i, 2], 0.5 * h0, tol,
                                     code, prm, d, poly, bu, bv, q3)
            if fj > f:
                f, n0, n1, n2 = fj, m0, m1, m2
        if mode == 1 or poly.size > 0:                         # polish: MCC creases (path) and the |n.d| crease of
                                                               # anisotropic strength can stall any local ascent;
            g, rad = 9, 0.5 * h0                               # re-seed from a 9 x 9 local grid while it improves
            for rep in range(4):
                ux, uy, uz, vx, vy, vz = _basis(n0, n1, n2)
                bf, b0, b1, b2 = f, n0, n1, n2
                for a in range(g):
                    for b in range(g):
                        sa_ = rad * (2.0 * a / (g - 1) - 1.0)
                        sb_ = rad * (2.0 * b / (g - 1) - 1.0)
                        c0 = n0 + sa_ * ux + sb_ * vx
                        c1 = n1 + sa_ * uy + sb_ * vy
                        c2 = n2 + sa_ * uz + sb_ * vz
                        r = np.sqrt(c0 * c0 + c1 * c1 + c2 * c2)
                        c0, c1, c2 = c0 / r, c1 / r, c2 / r
                        fc = _eval(mode, Ap, Bp, c0, c1, c2, code, prm, d, poly, bu, bv, q3)
                        if fc > bf:
                            bf, b0, b1, b2 = fc, c0, c1, c2
                bf, b0, b1, b2 = _refine(mode, Ap, Bp, bf, b0, b1, b2, rad / (g - 1), tol, code, prm, d, poly,
                                         bu, bv, q3)
                if bf <= f * (1.0 + 1e-12) + 1e-300:
                    break
                f, n0, n1, n2 = bf, b0, b1, b2
                rad *= 0.5
        if poly.size > 0:                                      # anisotropy: the maximum can sit ON the crease
            e1x, e1y, e1z, e2x, e2y, e2z = _basis(d[0], d[1], d[2])   # n.d = 0 (c = |n.d| has a kink): 1-D search on
            K = 360                                            # the great circle n = cos(phi) e1 + sin(phi) e2
            bphi, bf = 0.0, -np.inf
            for i in range(K):
                ph = np.pi * i / K
                fc = _eval(mode, Ap, Bp, np.cos(ph) * e1x + np.sin(ph) * e2x, np.cos(ph) * e1y + np.sin(ph) * e2y,
                           np.cos(ph) * e1z + np.sin(ph) * e2z, code, prm, d, poly, bu, bv, q3)
                if fc > bf:
                    bf, bphi = fc, ph
            w = np.pi / K
            while w > 0.25 * tol:                              # shrinking 1-D grid
                for j in (-2, -1, 1, 2):
                    ph = bphi + 0.5 * j * w
                    fc = _eval(mode, Ap, Bp, np.cos(ph) * e1x + np.sin(ph) * e2x, np.cos(ph) * e1y + np.sin(ph) * e2y,
                               np.cos(ph) * e1z + np.sin(ph) * e2z, code, prm, d, poly, bu, bv, q3)
                    if fc > bf:
                        bf, bphi = fc, ph
                w *= 0.5
            if bf > f:
                f, n0, n1, n2 = bf, np.cos(bphi) * e1x + np.sin(bphi) * e2x, np.cos(bphi) * e1y + np.sin(bphi) * e2y, \
                    np.cos(bphi) * e1z + np.sin(bphi) * e2z
        out_v[p] = _eval(mode, Ap, Bp, n0, n1, n2, code, prm, d, poly, bu, bv, q3)
        out_q[p, 0], out_q[p, 1], out_q[p, 2] = q3[0], q3[1], q3[2]
        sg = 1.0 if n2 >= 0 else -1.0
        out_n[p, 0], out_n[p, 1], out_n[p, 2] = sg * n0, sg * n1, sg * n2

    @numba.njit(parallel=True, cache=True)
    def _search_nb(mode, A, B, N, NB, q, h0, tol, code, prm, d, poly, out_v, out_q, out_n):
        P, M = A.shape[0], N.shape[0]
        T = A.shape[1] if mode == 1 else 1
        for p in numba.prange(P):
            _point(mode, A[p], B[p], N, NB, q, h0, tol, code, prm, d, poly, out_v, out_q, out_n, p, T, M)

def _run(mode, A, B, bound, crit, n_coarse, q, tol, screen, chunk):
    if numba is None:                                          # pragma: no cover
        raise ImportError("numba is required")
    P = len(A)
    val, qq, nn = np.full(P, np.nan), np.full((P, 3), np.nan), np.full((P, 3), np.nan)
    N, NB = coarse_grid(n_coarse)
    h0 = np.sqrt(2 * np.pi / n_coarse)
    order = np.argsort(-bound, kind="stable") if screen > 0 else np.arange(P)
    worst, done = -np.inf, 0
    step = chunk if screen > 0 else P
    for s in range(0, P, step):
        idx = order[s:s + step]
        if screen > 0 and bound[idx[0]] < screen * worst:
            break
        v, qb, nb = np.empty(len(idx)), np.empty((len(idx), 3)), np.empty((len(idx), 3))
        _search_nb(mode, np.ascontiguousarray(A[idx]), np.ascontiguousarray(B[idx]), N, NB, q, h0, tol,
                   crit.code, float(crit.prm), crit.d, crit.poly, v, qb, nb)
        val[idx], qq[idx], nn[idx] = v, qb, nb
        worst = max(worst, v.max())
        done += len(idx)
    return dict(value=val, tau_a=qq[:, 0], sn_a=qq[:, 1], sn_m=qq[:, 2], normal=nn, searched=~np.isnan(val),
                bound=bound, n_searched=done)


def _eig(S):
    return np.linalg.eigvalsh(to_matrix(S))                    # ascending


def search_prop(Sa, Sm, crit, n_coarse=800, q=4, tol=2e-4, screen=0.0, chunk=4096):
    """Constant-amplitude proportional loading.  Sa, Sm: (P, 6) amplitude and mean tensors.
    Returns dict of per-point arrays: value, tau_a, sn_a, sn_m, normal, searched, bound."""
    Sa = np.atleast_2d(np.asarray(Sa, float))
    Sm = np.atleast_2d(np.asarray(Sm, float))
    la, lm = _eig(Sa), _eig(Sm)
    TB = 0.5 * (la[:, 2] - la[:, 0])
    AB = np.abs(la).max(1)
    bound = crit.bound(TB, AB, lm[:, 2] + AB)
    return _run(0, Sa[:, None, :], Sm[:, None, :], bound, crit, n_coarse, q, tol, screen, chunk)


def search_path(H, crit, n_coarse=400, q=4, tol=2e-4, screen=0.0, chunk=512):
    """One cycle of a general (non-proportional) path.  H: (P, T, 6)."""
    H = np.asarray(H, float)
    if H.ndim == 2:
        H = H[None]
    C = 0.5 * (H.max(1) + H.min(1))
    l = _eig(H - C[:, None, :])                                # (P, T, 3)
    TB = 0.5 * (l[..., 2] - l[..., 0]).max(1)
    AB = np.abs(l).max((1, 2))
    NB = _eig(H)[..., 2].max(1)
    bound = crit.bound(TB, AB, NB)
    return _run(1, H, np.zeros((len(H), 1, 6)), bound, crit, n_coarse, q, tol, screen, chunk)


def superpose(unit_fields, loads):
    """Linear superposition: unit_fields (K, P, 6) stresses of K unit load cases, loads (T, K) -> history (P, T, 6)."""
    return np.einsum("tk,kpj->ptj", np.asarray(loads, float), np.asarray(unit_fields, float))
