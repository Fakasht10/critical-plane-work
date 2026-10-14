"""Speed / accuracy benchmark of the critical-plane kernel (cplane.py).

    python cplane_bench.py   -> results/CPLANE_SPEED.md, cplane_speed.json, cplane_speed.png

Methods (Findley k = 0.3, constant-amplitude proportional and a 32-sample non-proportional cycle):
  naive   : Python loop over points and a 2 deg (theta, phi) grid (the usual script pattern), timed on a subset and
            scaled linearly;
  dense   : vectorised GEMM on a 1 deg hemisphere grid (20 600 normals), no refinement;
  ours    : defaults, 800 (prop) / 400 (path) coarse normals, 4 best coarse local maxima refined to 2e-4 rad;
  ours+scr: same with hot-spot screening at 0.5 x worst damage.
Reference for errors: prop - ours with 8000 coarse normals, 8 candidates, tol 1e-7 (itself checked against a 60 000-
normal brute force in tests/test_cplane.py); path - independent brute force on 200 000 normals (first 200 points).
Synthetic fields: random tensors scaled by a log-normal factor (sigma 1.0) so that, as in a notched FE model, few
points carry high stress.  The FE-field timing is in CCX_NGSOLVE.md.
"""
import json
import os
import time

import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

import cplane as cp
from plots import _save

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "results")
K = 0.3


def field(P, seed=0):
    rng = np.random.default_rng(seed)
    s = rng.lognormal(0, 1.0, (P, 1))
    return rng.normal(size=(P, 6)) * s, 0.3 * rng.normal(size=(P, 6)) * s


def path(P, T=32, seed=0):
    Sa, Sm = field(P, seed)
    rng = np.random.default_rng(seed + 1)
    Sb = rng.normal(size=(P, 6)) * np.abs(Sa).max(1, keepdims=True) * 0.5
    t = np.linspace(0, 2 * np.pi, T, endpoint=False)
    return Sm[:, None] + np.sin(t)[None, :, None] * Sa[:, None] + np.cos(t)[None, :, None] * Sb[:, None]


def naive_prop(Sa, Sm, step_deg=2.0):
    th = np.radians(np.arange(0, 90 + 1e-9, step_deg))
    ph = np.radians(np.arange(0, 360, step_deg))
    out = np.empty(len(Sa))
    for p in range(len(Sa)):
        A, M = cp.to_matrix(Sa[p]), cp.to_matrix(Sm[p])
        best = -np.inf
        for t in th:
            for f in ph:
                n = np.array([np.sin(t) * np.cos(f), np.sin(t) * np.sin(f), np.cos(t)])
                tr, trm = A @ n, M @ n
                sn = n @ tr
                tau = np.linalg.norm(tr - sn * n)
                best = max(best, tau + K * (abs(sn) + n @ trm))
        out[p] = best
    return out


def dense(fn, data, m=20600):
    N = cp.fib_hemisphere(m)
    W = cp.proj_rows(N)
    out = []
    for s in range(0, len(data[0]), 256 if fn is cp._quant_path else 2048):
        ta, sa, sm = fn(W, *(x[s:s + 256 if fn is cp._quant_path else s + 2048] for x in data))
        out.append((ta + K * (sa + sm)).max(1))
    return np.concatenate(out)


_BRUTE = {}


def brute_path(H, m=200_000):
    key = H.tobytes()
    if key not in _BRUTE:
        N = cp.fib_hemisphere(m)
        W = cp.proj_rows(N)
        _BRUTE[key] = np.array([(lambda q: (q[0] + K * (q[1] + q[2])).max())(cp._quant_path(W, H[p:p + 1]))
                                for p in range(len(H))])
    return _BRUTE[key]


def timed(f, *a, **kw):
    t = time.perf_counter()
    r = f(*a, **kw)
    return r, time.perf_counter() - t


def rel(a, ref):
    m = ~np.isnan(a)
    return float(np.max(np.abs(a[m] - ref[m]) / np.abs(ref[m])))


def main():
    crit = cp.Findley(K)
    cp.search_prop(*field(10), crit)
    cp.search_path(path(10), crit)                              # JIT warm-up
    rows = []
    for P in (1_000, 10_000, 100_000, 1_000_000):
        Sa, Sm = field(P)
        sub = min(P, 20_000)
        ref, _ = timed(cp.search_prop, Sa[:sub], Sm[:sub], crit, n_coarse=8000, q=8, tol=1e-7)
        ref = ref["value"]
        r, t_ours = timed(cp.search_prop, Sa, Sm, crit)
        rs, t_scr = timed(cp.search_prop, Sa, Sm, crit, screen=0.5)
        row = dict(case="prop", P=P, ours_s=t_ours, ours_err=rel(r["value"][:sub], ref), scr_s=t_scr,
                   scr_searched=int(rs["n_searched"]), scr_max_ok=bool(np.isclose(np.nanmax(rs["value"]), r["value"].max())))
        nn = 30
        _, t = timed(naive_prop, Sa[:nn], Sm[:nn])
        row["naive_s"] = t * P / nn
        row["naive_err"] = rel(naive_prop(Sa[:nn], Sm[:nn]), ref[:nn])
        if P <= 100_000:
            d, row["dense_s"] = timed(dense, cp._quant_prop, (Sa, Sm))
            row["dense_err"] = rel(d[:sub], ref)
        rows.append(row)
        print(row, flush=True)
    for P in (1_000, 10_000, 100_000):
        H = path(P)
        sub = 200
        ref = brute_path(H[:sub])
        r, t_ours = timed(cp.search_path, H, crit)
        rs, t_scr = timed(cp.search_path, H, crit, screen=0.5)
        row = dict(case="path T=32", P=P, ours_s=t_ours, ours_err=rel(r["value"][:sub], ref), scr_s=t_scr,
                   scr_searched=int(rs["n_searched"]), scr_max_ok=bool(np.isclose(np.nanmax(rs["value"]), r["value"].max())))
        if P <= 1_000:
            d, row["dense_s"] = timed(dense, cp._quant_path, (H,))
            row["dense_err"] = rel(d[:sub], ref)
        rows.append(row)
        print(row, flush=True)
    json.dump(rows, open(os.path.join(OUT, "cplane_speed.json"), "w"), indent=1)
    return rows


def report(rows=None):
    rows = rows or json.load(open(os.path.join(OUT, "cplane_speed.json")))
    f = lambda x: "-" if x is None else ("%.3g" % x)
    L = ["# Critical-plane kernel: speed and accuracy (cplane.py)", "",
         "Findley (k = 0.3). Synthetic log-normal-scaled stress fields; 4 CPU cores, FP64. Times in seconds, errors are "
         "the max relative error of the damage parameter vs a converged reference (prop: 8000 coarse normals, 8 candidates, "
         "tol 1e-7; path: brute force on 200 000 normals, first 200 points). Defaults: 800 (prop) / 400 (path) coarse "
         "normals, 4 local maxima refined to 2e-4 rad. `naive` = Python loop, 2 deg grid, timed on 30 points and scaled. `dense` = vectorised 1 deg grid "
         "(20 600 normals). `screen 0.5` = exact search only where the upper bound exceeds 0.5 x worst damage.", "",
         "| loading | points | naive [s] | naive err | dense 1deg [s] | dense err | **ours [s]** | **ours err** | "
         "ours+screen [s] | points searched | worst point kept |", "|" + "---|" * 11]
    for r in rows:
        L.append("| %s | %s | %s | %s | %s | %s | **%s** | **%s** | %s | %s | %s |" % (
            r["case"], f"{r['P']:,}", f(r.get("naive_s")), f(r.get("naive_err")), f(r.get("dense_s")), f(r.get("dense_err")),
            f(r["ours_s"]), f(r["ours_err"]), f(r["scr_s"]), f"{r['scr_searched']:,}", "yes" if r["scr_max_ok"] else "NO"))
    pr = [r for r in rows if r["case"] == "prop"]
    L += ["", "Throughput (proportional, full field): %.2g points/s; naive: %.2g points/s." % (
        pr[-1]["P"] / pr[-1]["ours_s"], pr[-1]["P"] / pr[-1]["naive_s"]), "",
        "![speed](cplane_speed.png)"]
    open(os.path.join(OUT, "CPLANE_SPEED.md"), "w").write("\n".join(L) + "\n")
    fig, ax = plt.subplots(1, 2, figsize=(9, 3.4))
    for case, mk in (("prop", "o"), ("path T=32", "s")):
        rr = [r for r in rows if r["case"] == case]
        P = [r["P"] for r in rr]
        ax[0].loglog(P, [r["ours_s"] for r in rr], mk + "-", label=f"ours, {case}")
        ax[0].loglog(P, [r["scr_s"] for r in rr], mk + "--", label=f"ours + screen, {case}")
        d = [(r["P"], r["dense_s"]) for r in rr if "dense_s" in r]
        ax[0].loglog(*zip(*d), mk + ":", label=f"dense 1 deg, {case}")
        if case == "prop":
            ax[0].loglog(P, [r["naive_s"] for r in rr], mk + "-.", color="k", label="naive loop (scaled)")
        ax[1].semilogx(P, [r["ours_err"] for r in rr], mk + "-", label=f"ours, {case}")
        ax[1].semilogx(*zip(*[(r["P"], r["dense_err"]) for r in rr if "dense_err" in r]), mk + ":", label=f"dense, {case}")
    ax[0].set(xlabel="points", ylabel="wall time [s]", title="(a) time")
    ax[1].set(xlabel="points", ylabel="max rel. error", yscale="log", title="(b) accuracy")
    for a in ax:
        a.grid(alpha=0.3, which="both")
    ax[0].legend(fontsize=6.5, loc="lower right")
    ax[1].legend(fontsize=7, loc="center right")
    fig.tight_layout()
    _save(fig, OUT, "cplane_speed")


if __name__ == "__main__":
    report(main())
