"""Verification of the critical-plane kernel: closed forms, MCC, brute force, invariance, screening."""
import numpy as np
import pytest

import cplane as cp


def v(xx=0, yy=0, zz=0, xy=0, yz=0, zx=0):
    return np.array([[xx, yy, zz, xy, yz, zx]], float)


def rand_rot(rng):
    Q, R = np.linalg.qr(rng.normal(size=(3, 3)))
    return Q * np.sign(np.diag(R))


def rotate(S, Q):
    return cp.to_voigt(Q @ cp.to_matrix(S) @ Q.T)


@pytest.mark.parametrize("k", [0.0, 0.2, 0.5])
def test_findley_uniaxial_and_torsion_closed_form(k):
    sa = 10.0
    r = cp.search_prop(v(xx=sa), v(), cp.Findley(k))
    assert r["value"][0] == pytest.approx(sa / 2 * (k + np.sqrt(1 + k * k)), rel=1e-8)
    r = cp.search_prop(v(xy=sa), v(), cp.Findley(k))
    assert r["value"][0] == pytest.approx(sa * np.sqrt(1 + k * k), rel=1e-8)


def test_mwcm_uniaxial_plane_and_rho():
    r = cp.search_prop(v(xx=10.0), v(xx=10.0), cp.MWCM())          # R = 0
    assert r["tau_a"][0] == pytest.approx(5.0, rel=1e-9)
    rho = (r["sn_m"] + r["sn_a"]) / r["tau_a"]
    assert rho[0] == pytest.approx(2.0, rel=1e-4)                   # sn_max = 10 on the 45 deg plane
    assert abs(r["normal"][0, 0]) == pytest.approx(np.sqrt(0.5), abs=2e-3)


def test_mcc_circle_ellipse_collinear():
    rng = np.random.default_rng(1)
    t = rng.permutation(np.linspace(0, 2 * np.pi, 97, endpoint=False))
    assert cp.mcc_radius(3 * np.cos(t)[None], 3 * np.sin(t)[None] + 1)[0] == pytest.approx(3.0, rel=1e-12)
    assert cp.mcc_radius(5 * np.cos(t)[None], 2 * np.sin(t)[None])[0] == pytest.approx(5.0, rel=1e-3)
    x = rng.uniform(-1, 4, 50)
    assert cp.mcc_radius(x[None], (2 * x)[None])[0] == pytest.approx(np.ptp(x) * np.sqrt(5) / 2, rel=1e-12)


def test_mcc_against_minimax():
    from scipy.optimize import minimize
    rng = np.random.default_rng(2)
    for _ in range(20):
        P = rng.normal(size=(rng.integers(3, 40), 2)) * rng.uniform(0.1, 10)
        r = cp.mcc_radius(P[None, :, 0], P[None, :, 1])[0]
        f = lambda c: np.max(np.hypot(*(P - c).T))
        ref = minimize(f, P.mean(0), method="Nelder-Mead", options=dict(xatol=1e-12, fatol=1e-14, maxiter=20000)).fun
        assert r == pytest.approx(ref, rel=1e-6)
        assert np.all(np.hypot(*(P - P.mean(0)).T) <= np.inf)


def test_out_of_phase_tension_torsion():
    """sigma = s sin t, tau = s/2 cos t (90 deg): on the 45-deg-family planes the analytic Findley optimum is checked
    against a dense brute-force scan of the same kernel quantities."""
    s, k = 10.0, 0.3
    t = np.linspace(0, 2 * np.pi, 64, endpoint=False)
    H = np.zeros((1, t.size, 6))
    H[0, :, 0] = s * np.sin(t)
    H[0, :, 3] = s / 2 * np.cos(t)
    r = cp.search_path(H, cp.Findley(k))
    N = cp.fib_hemisphere(40000)
    ta, sa, sm = cp._quant_path(cp.proj_rows(N), H)
    brute = np.max(ta + k * (sa + sm))
    assert r["value"][0] >= brute - 1e-9
    assert r["value"][0] == pytest.approx(brute, rel=2e-4)


def test_path_equals_prop_for_proportional_history():
    rng = np.random.default_rng(3)
    Sa, Sm = rng.normal(size=(50, 6)), rng.normal(size=(50, 6))
    t = np.linspace(0, 2 * np.pi, 33)
    H = Sm[:, None, :] + np.sin(t)[None, :, None] * Sa[:, None, :]
    crit = cp.Findley(0.25)
    a, b = cp.search_prop(Sa, Sm, crit), cp.search_path(H, crit)
    np.testing.assert_allclose(a["value"], b["value"], rtol=1e-6)


@pytest.mark.parametrize("crit", [cp.Findley(0.3), cp.SWTStress(), cp.MWCM()])
def test_against_dense_brute_force(crit):
    rng = np.random.default_rng(4)
    Sa, Sm = rng.normal(size=(200, 6)), 0.5 * rng.normal(size=(200, 6))
    r = cp.search_prop(Sa, Sm, crit)
    N = cp.fib_hemisphere(60000)
    ta, sa, sm = cp._quant_prop(cp.proj_rows(N), Sa, Sm)
    brute = crit.objective(ta, sa, sm, N).max(1)
    assert np.all(r["value"] >= brute - 1e-9)                       # refinement never worse than a 0.6 deg grid
    assert np.max((r["value"] - brute) / np.abs(brute)) < 2e-4


def test_rotation_invariance_with_build_direction():
    rng = np.random.default_rng(5)
    Sa, Sm = rng.normal(size=(40, 6)), rng.normal(size=(40, 6))
    Q = rand_rot(rng)
    weak = (1.0, 0.0, -0.4)                                         # inter-layer planes (n || d) 40 % weaker
    d = np.array([0.0, 0.0, 1.0])
    a = cp.search_prop(Sa, Sm, cp.Findley(0.3, build_dir=d, strength_poly=weak))
    b = cp.search_prop(rotate(Sa, Q), rotate(Sm, Q), cp.Findley(0.3, build_dir=Q @ d, strength_poly=weak))
    np.testing.assert_allclose(a["value"], b["value"], rtol=1e-6)            # convergence level (tol 2e-4 rad)
    iso = cp.search_prop(Sa, Sm, cp.Findley(0.3))
    assert np.all(a["value"] >= iso["value"] - 1e-12)               # weaker planes can only raise damage


def test_bound_is_upper_bound_and_screening_finds_hotspot():
    rng = np.random.default_rng(6)
    P = 3000
    Sa = rng.normal(size=(P, 6)) * rng.lognormal(0, 0.6, (P, 1))
    Sm = 0.3 * rng.normal(size=(P, 6))
    crit = cp.Findley(0.3, build_dir=[0, 0, 1], strength_poly=(1.0, -0.3))
    full = cp.search_prop(Sa, Sm, crit, chunk=500)
    assert np.all(full["bound"] >= full["value"] - 1e-9)
    scr = cp.search_prop(Sa, Sm, crit, screen=0.7, chunk=500)
    assert np.nanmax(scr["value"]) == pytest.approx(full["value"].max(), rel=1e-12)
    keep = full["value"] >= 0.7 * full["value"].max()
    assert np.all(scr["searched"][keep])                            # every point within 70 % of the worst is exact
    assert scr["n_searched"] < P

    H = Sm[:200, None, :] + np.sin(np.linspace(0, 6.28, 16))[None, :, None] * Sa[:200, None, :]
    H[:, :, 3] += np.cos(np.linspace(0, 6.28, 16))[None, :] * Sa[:200, None, 0]
    r = cp.search_path(H, crit)
    assert np.all(r["bound"] >= r["value"] - 1e-9)


def test_superpose_and_abaqus_order():
    rng = np.random.default_rng(7)
    U = rng.normal(size=(2, 5, 6))
    L = rng.normal(size=(9, 2))
    H = cp.superpose(U, L)
    np.testing.assert_allclose(H[3, 4], L[4, 0] * U[0, 3] + L[4, 1] * U[1, 3])
    a = np.array([1, 2, 3, 12, 13, 23.0])
    np.testing.assert_allclose(cp.from_abaqus(a), [1, 2, 3, 12, 23, 13])


def test_walker_limits_and_closed_form():
    rng = np.random.default_rng(8)
    Sa, Sm = rng.normal(size=(60, 6)), rng.normal(size=(60, 6))
    np.testing.assert_allclose(cp.search_prop(Sa, Sm, cp.WalkerStress(0.5))["value"],
                               cp.search_prop(Sa, Sm, cp.SWTStress())["value"], rtol=1e-6)
    smax = np.maximum(np.linalg.eigvalsh(cp.to_matrix(Sm + Sa))[:, 2],          # g = 0: max principal over the two
                      np.linalg.eigvalsh(cp.to_matrix(Sm - Sa))[:, 2])          # peak tensors of the cycle
    np.testing.assert_allclose(cp.search_prop(Sa, Sm, cp.WalkerStress(0.0))["value"], smax, rtol=1e-6)
    for g in (0.0, 0.3, 1.0):                                                   # uniaxial: R = -1 -> s_a; R = 0 -> 2^(1-g) s_a
        assert cp.search_prop(v(xx=4.0), v(), cp.WalkerStress(g))["value"][0] == pytest.approx(4.0, rel=1e-8)
        assert cp.search_prop(v(xx=4.0), v(xx=4.0), cp.WalkerStress(g))["value"][0] == pytest.approx(4.0 * 2 ** (1 - g), rel=1e-8)
    r = cp.search_prop(Sa, Sm, cp.WalkerStress(0.3))
    assert np.all(r["bound"] >= r["value"] - 1e-9)


def test_anisotropic_crease_prop_against_brute_force():
    """strength poly with a positive linear term puts a crease on n.d = 0 (audit finding): polish must reach it."""
    rng = np.random.default_rng(11)
    Sa, Sm = rng.normal(size=(150, 6)), 0.4 * rng.normal(size=(150, 6))
    crit = cp.Findley(0.3, build_dir=[0.3, -0.5, 0.8], strength_poly=(1.0, 0.8, -1.2, 0.6))
    r = cp.search_prop(Sa, Sm, crit)
    N = cp.fib_hemisphere(200000)
    d = crit.d
    e1 = np.cross(d, [1.0, 0, 0])
    e1 /= np.linalg.norm(e1)
    ph = np.linspace(0, np.pi, 100001)                                  # the crease circle n . d = 0, densely
    N = np.vstack([N, np.cos(ph)[:, None] * e1 + np.sin(ph)[:, None] * np.cross(d, e1)])
    ta, sa, sm = cp._quant_prop(cp.proj_rows(N), Sa, Sm)
    brute = crit.objective(ta, sa, sm, N).max(1)
    assert np.all(r["value"] >= brute * (1 - 1e-8))
    assert np.max((r["value"] - brute) / np.abs(brute)) < 1e-4


def test_exact_strength_minimum_keeps_bound_valid():
    crit = cp.WalkerStress(0.7, build_dir=[0, 0, 1], strength_poly=(2.0, -3.0, 1.6))
    assert crit.min_strength() <= 0.59375
    rng = np.random.default_rng(12)
    Sa, Sm = rng.normal(size=(2000, 6)), rng.normal(size=(2000, 6))
    r = cp.search_prop(Sa, Sm, crit)
    assert np.all(r["bound"] >= r["value"])


def test_mcc_large_offset():
    rng = np.random.default_rng(13)
    t = rng.permutation(np.linspace(0, 2 * np.pi, 50, endpoint=False))
    r = cp.mcc_radius((1e8 + 1e-3 * np.cos(t))[None], (-3e7 + 1e-3 * np.sin(t))[None])[0]
    assert r == pytest.approx(1e-3, rel=5e-5)      # inputs at 1e8 are quantised to ~1.5e-8 (1.5e-5 of the radius)
