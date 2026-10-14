"""Non-proportional multiaxial fatigue of CFRP tubes (Duda et al., ACME 25:212, 2025): critical-plane correlation
and leave-one-set-out prediction with cplane.

    python cp_duda_cfrp.py   -> results/CP_DUDA_CFRP.md, cp_duda_cfrp.json, cp_duda_cfrp.png

Data: duda_data.CFRP (Table 4, 104 tests, 9 sets = R {-1, 0.05, 0.5} x phase {0, 45, 90 deg}, tau_xy,a = sigma_x,a).

Stress analysis (transparent, linear): thin-walled [+-30] tube, classical laminate theory (membrane), lamina
properties of Table 1; ply stresses (sigma_1, sigma_2, tau_12) histories over one cycle (64 samples) from
sigma_x(t) = sigma_m + sigma_a sin t, tau_xy(t) = tau_m + tau_a sin(t - delta), same R on both channels.
The paper used a nonlinear FE critical point instead (zig-zag area); a constant stress-concentration factor
does not change correlations, the nonlinear shear response does (their lambda_12 depends on R; CLT's does not).

Damage parameters (per specimen; max over the two plies):
  PAPER     d_tau12 + k sigma_2,max on the fibre-parallel plane, k = 0.2 (paper's value);
  CP-FIBRE  cplane.search_path, Findley tau_a + k sigma_n,max, planes restricted to contain the fibre (strength
            polynomial 1 + 1e4 c^2, c = |n . e1|; sigma_1 dropped: it has no traction on such planes). Shear
            amplitude = minimum circumscribed circle. With k = 0.1 it must equal PAPER / 2 (check);
            k fitted (max r) on the positive-R sets;
  CP-ISO    cplane.search_path on the applied (geometric) stresses, isotropic Findley, k fitted the same way:
            what a metals-style analysis would do.
Statistics: Basquin log10 N = a + b log10 P on failures; Pearson r, rms log10 N; leave-one-set-out (fit on 8 sets,
predict the 9th) rms log10 N per set.
"""
import json
import os

import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

import cplane as cp
import duda_data as D
from plots import _save

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "results")
T_SAMPLES = 64
POS = ("C1", "C2", "C4", "C5", "C7", "C8")                               # positive R
PAPER_FIT = D.PAPER["fit_sets"]


# ------------------------------------------------------------------------------------------------------- CLT
def q_matrices():
    p = D.PLY
    nu21 = p["nu12"] * p["E2"] / p["E1"]
    d = 1 - p["nu12"] * nu21
    Q = np.array([[p["E1"] / d, p["nu12"] * p["E2"] / d, 0], [p["nu12"] * p["E2"] / d, p["E2"] / d, 0], [0, 0, p["G12"]]])

    def T_sig(th):                                  # global (x, y, xy) stress -> material (1, 2, 12)
        c, s = np.cos(th), np.sin(th)
        return np.array([[c * c, s * s, 2 * c * s], [s * s, c * c, -2 * c * s], [-c * s, c * s, c * c - s * s]])

    def Qbar(th):
        Ts = T_sig(th)
        R = np.diag([1, 1, 2.0])
        Te = R @ Ts @ np.linalg.inv(R)              # engineering-strain transform
        return np.linalg.inv(Ts) @ Q @ Te
    ths = np.radians([D.PLY["angle"], -D.PLY["angle"]])
    A = sum(Qbar(t) for t in ths) / len(ths)
    return ths, Qbar, T_sig, A


THS, QBAR, TSIG, ABAR = q_matrices()
COMP = np.linalg.inv(ABAR)


def ply_stresses(sx, txy):
    """Membrane stresses (sx, 0, txy) [..] -> list over plies of material stresses (s1, s2, t12) [.., 3]."""
    sg = np.stack([sx, np.zeros_like(sx), txy], -1)
    eps = sg @ COMP.T
    return [(eps @ QBAR(th).T) @ TSIG(th).T for th in THS]


def history(sa, R, delta):
    t = np.linspace(0, 2 * np.pi, T_SAMPLES, endpoint=False)
    sm = sa * (1 + R) / (1 - R)
    return sm + sa * np.sin(t), sm + sa * np.sin(t - np.radians(delta))   # tau_a = sigma_a (lambda_T = 1)


# ------------------------------------------------------------------------------------------------ parameters
def specimens():
    rows = []
    for s, (delta, R, data) in D.CFRP.items():
        for sa, N in data:
            rows.append(dict(set=s, delta=delta, R=R, sa=sa, N=N, runout=N >= D.RUNOUT))
    return rows


def param_paper(rows, k=D.PAPER["k"]):
    out = []
    for r in rows:
        sx, txy = history(r["sa"], r["R"], r["delta"])
        out.append(max(np.ptp(p[:, 2]) + k * p[:, 1].max() for p in ply_stresses(sx, txy)))
    return np.array(out)


def ply_tensor_histories(rows):
    """[n, T, 6] per ply, material axes, sigma_1 dropped (no traction on fibre-containing planes)."""
    H = {0: [], 1: []}
    for r in rows:
        sx, txy = history(r["sa"], r["R"], r["delta"])
        for i, p in enumerate(ply_stresses(sx, txy)):
            h = np.zeros((T_SAMPLES, 6))
            h[:, 1], h[:, 3] = p[:, 1], p[:, 2]                                 # yy = sigma_2, xy = tau_12
            H[i].append(h)
    return {i: np.array(v) for i, v in H.items()}


def geometric_histories(rows):
    H = np.zeros((len(rows), T_SAMPLES, 6))
    for j, r in enumerate(rows):
        sx, txy = history(r["sa"], r["R"], r["delta"])
        H[j, :, 0], H[j, :, 3] = sx, txy
    return H


def param_cp_fibre(Hply, k):
    crit = cp.Findley(k, build_dir=[1, 0, 0], strength_poly=(1.0, 0.0, 1e4))
    return np.maximum(*(cp.search_path(Hply[i], crit)["value"] for i in (0, 1)))


def param_cp_iso(Hg, k):
    return cp.search_path(Hg, cp.Findley(k))["value"]


# ------------------------------------------------------------------------------------------------ statistics
def fit(P, N, mask):
    b, a = np.polyfit(np.log10(P[mask]), np.log10(N[mask]), 1)
    res = np.log10(N[mask]) - (a + b * np.log10(P[mask]))
    r = np.corrcoef(np.log10(P[mask]), np.log10(N[mask]))[0, 1]
    return dict(a=float(a), b=float(b), r=float(r), rms=float(np.sqrt(np.mean(res ** 2))))


def loso(P, N, sets, fail):
    out = {}
    for s in sorted(set(sets)):
        tr, te = fail & (sets != s), fail & (sets == s)
        f = fit(P, N, tr)
        res = np.log10(N[te]) - (f["a"] + f["b"] * np.log10(P[te]))
        out[s] = dict(rms=float(np.sqrt(np.mean(res ** 2))), bias=float(res.mean()))
    return out


def loso_nested(Pk, N, sets, fail, pos_sets):
    """Leave-one-set-out with the constant k re-selected inside each fold (on the training positive-R sets only)."""
    out = {}
    for s in sorted(set(sets)):
        tr = fail & (sets != s)
        sel = tr & np.isin(sets, pos_sets)
        k = min(Pk, key=lambda kk: fit(Pk[kk], N, sel)["r"])
        f = fit(Pk[k], N, tr)
        te = fail & (sets == s)
        res = np.log10(N[te]) - (f["a"] + f["b"] * np.log10(Pk[k][te]))
        out[s] = dict(rms=float(np.sqrt(np.mean(res ** 2))), bias=float(res.mean()), k=float(k))
    return out


def oracle_r(sets_sel):
    """Best |r| any load-proportional parameter can reach: one free scale factor per set (a homogeneous parameter at
    fixed R and phase is c_set x sigma_a), optimised over the scales."""
    from scipy.optimize import minimize
    data = [(i, np.log10(a), np.log10(N)) for i, s in enumerate(sets_sel) for a, N in D.CFRP[s][2] if N < D.RUNOUT]
    i, x, y = (np.array(v) for v in zip(*data))
    f = lambda c: np.corrcoef(x + np.r_[0, c][i.astype(int)], y)[0, 1]
    r = minimize(f, np.zeros(len(sets_sel) - 1), method="Nelder-Mead", options=dict(maxiter=40000, xatol=1e-9, fatol=1e-12))
    return float(-r.fun), (10 ** np.r_[0, r.x]).tolist()


def within_set_r():
    out = {}
    for s, (delta, R, data) in D.CFRP.items():
        d = np.array([(a, N) for a, N in data if N < D.RUNOUT])
        out[s] = float(np.corrcoef(np.log10(d[:, 0]), np.log10(d[:, 1]))[0, 1])
    return out


def walker_scan(Hg, N, mask):
    return [(float(g), float(np.corrcoef(np.log10(cp.search_path(Hg, cp.WalkerStress(g))["value"][mask]),
                                          np.log10(N[mask]))[0, 1])) for g in (0.0, 0.25, 0.5, 0.75, 1.0)]


def main():
    rows = specimens()
    N = np.array([r["N"] for r in rows], float)
    sets = np.array([r["set"] for r in rows])
    fail = ~np.array([r["runout"] for r in rows])
    pos = np.isin(sets, POS) & fail
    pfit = np.isin(sets, PAPER_FIT) & fail
    Hply, Hg = ply_tensor_histories(rows), geometric_histories(rows)
    out = dict(sets={s: dict(delta=D.CFRP[s][0], R=D.CFRP[s][1]) for s in D.CFRP})
    # linear CLT lambda_12 vs the paper's nonlinear FE values
    p1 = ply_stresses(np.array([1.0]), np.array([1.0]))                  # unit sigma_x = tau_xy (linear: R-independent)
    out["ply_unit"] = {f"{d:+.0f}": dict(s1=float(q[0, 0]), s2=float(q[0, 1]), t12=float(q[0, 2]),
                                         lam=float(abs(q[0, 2] / q[0, 1])))
                       for d, q in zip(np.degrees(THS), p1)}
    P = {}
    P["nominal"] = np.array([r["sa"] for r in rows])
    P["paper"] = param_paper(rows)
    chk = param_cp_fibre(Hply, 0.1)
    out["check_cp_vs_paper"] = float(np.max(np.abs(chk / (P["paper"] / 2) - 1)))
    print("cplane fibre-plane (k=0.1) vs paper/2: max rel diff %.2e" % out["check_cp_vs_paper"], flush=True)
    ks = np.round(np.arange(0.0, 1.51, 0.05), 3)
    Pk_f = {float(k): param_cp_fibre(Hply, k) for k in ks}
    Pk_i = {float(k): param_cp_iso(Hg, k) for k in ks}
    scan_f = [(k, fit(v, N, pos)["r"]) for k, v in Pk_f.items()]
    kf = min(scan_f, key=lambda x: x[1])[0]                                 # r is negative: most negative = best
    scan_i = [(k, fit(v, N, pos)["r"]) for k, v in Pk_i.items()]
    ki = min(scan_i, key=lambda x: x[1])[0]
    P["cp_fibre"] = param_cp_fibre(Hply, kf)
    P["cp_iso"] = param_cp_iso(Hg, ki)
    out["k"] = dict(cp_fibre=kf, cp_iso=ki, scan_fibre=scan_f, scan_iso=scan_i)
    res = {}
    for name, p in P.items():
        res[name] = dict(paper_subset=fit(p, N, pfit), positive_R=fit(p, N, pos), all=fit(p, N, fail),
                         loso=loso(p, N, sets, fail), P=p.tolist())
        print(name, {kk: (round(v["r"], 3), round(v["rms"], 3)) for kk, v in res[name].items() if kk in
                     ("paper_subset", "positive_R", "all")}, flush=True)
    res["cp_fibre"]["loso"] = loso_nested(Pk_f, N, sets, fail, POS)        # k re-selected per fold (no leakage)
    res["cp_iso"]["loso"] = loso_nested(Pk_i, N, sets, fail, POS)
    out["res"] = res
    out["oracle"] = dict(paper_subset=oracle_r(PAPER_FIT), positive_R=oracle_r(POS), all=oracle_r(tuple(D.CFRP)))
    out["within_r"] = within_set_r()
    out["walker_iso_scan_all"] = walker_scan(Hg, N, fail)
    out["N"], out["setcol"], out["fail"] = N.tolist(), sets.tolist(), fail.tolist()
    json.dump(out, open(os.path.join(OUT, "cp_duda_cfrp.json"), "w"), indent=1)
    return out


LAB = {"nominal": "nominal σx,a", "paper": "paper FS (CLT stresses, k = 0.2)",
       "cp_fibre": "cplane, fibre planes (Findley)", "cp_iso": "cplane, isotropic (Findley)"}


def report(out=None):
    out = out or json.load(open(os.path.join(OUT, "cp_duda_cfrp.json")))
    res = out["res"]
    L = ["# Non-proportional multiaxial fatigue of CFRP tubes (Duda et al. 2025): critical-plane correlation", "",
         "Data: Duda, Smolnicki, Lesiuk, Fatemi, Arch. Civ. Mech. Eng. 25 (2025) 212 (CC BY), Table 4: 104 thin-walled "
         "[+-30] CFRP tubes, 9 sets = R (-1, 0.05, 0.5) x phase (0, 45, 90 deg), tau_xy,a = sigma_x,a, 14 run-outs. "
         "Ply stresses by linear CLT (Table 1 lamina properties); the paper used a nonlinear FE critical point.", "",
         "Check: cplane on fibre-containing planes with Findley k = 0.1 reproduces the paper's parameter / 2 to "
         "%.1e (the shear path on the fibre plane is 1D, so MCC = range / 2)." % out["check_cp_vs_paper"],
         "Linear CLT, unit sigma_x = tau_xy (independent of R): " + "; ".join(
             "ply %s: sigma_1 %.2f, sigma_2 %.3f, tau_12 %.3f (|tau_12/sigma_2| %.2f)" % (k, v["s1"], v["s2"], v["t12"], v["lam"])
             for k, v in out["ply_unit"].items()) + ". **sigma_2 is compressive in both plies** (laminate Poisson "
         "effect), so k sigma_2,max <= 0 for positive R and the paper's parameter cannot work with CLT stresses; the "
         "fibre-plane fit of k hits its lower bound 0 (Findley requires k >= 0) - a constrained, not a true, optimum. "
         "The paper's nonlinear FE critical point has lambda_12 = 0.35 / 0.25 / 0.18 (R = -1 / 0.05 / 0.5) with tensile "
         "sigma_2. Membrane-only CLT also ignores the B16/B26 coupling of the two-layer +-30 wall and the ~25 % "
         "through-wall variation of torsional shear (t / R_m = 0.22).", "",
         "Fitted Findley k (max |r| on positive-R sets, in-sample): fibre planes %.2f, isotropic %.2f." % (
             out["k"]["cp_fibre"], out["k"]["cp_iso"]), "",
         "## Correlation (Basquin on failures): Pearson r / rms log10 N (in-sample for the fitted-k rows)", "",
         "| parameter | paper subset C1 C2 C7 C8 | positive R (6 sets) | all 9 sets |", "|---|---|---|---|"]
    for n, q in res.items():
        L.append("| %s | %.2f / %.2f | %.2f / %.2f | %.2f / %.2f |" % (
            LAB[n], q["paper_subset"]["r"], q["paper_subset"]["rms"], q["positive_R"]["r"], q["positive_R"]["rms"],
            q["all"]["r"], q["all"]["rms"]))
    o = out["oracle"]
    L.append("| *ceiling: best \\|r\\| of ANY load-proportional parameter (free scale per set)* | *%.2f* | *%.2f* | *%.2f* |" % (
        o["paper_subset"][0], o["positive_R"][0], o["all"][0]))
    t3 = D.PAPER["table3_r"]
    L += ["", "Paper (nonlinear FE stresses, C1 C2 C7 C8): r = 0.84.", "",
          "**Why every correlation is weak.** Within a set the stress range is narrow (e.g. 24-30 MPa) and life scatters "
          "over 1-2 decades. Any parameter proportional to the load at fixed R and phase (all stress-based critical-plane "
          "parameters with linear stresses: Findley, MCC, Walker are degree-1 homogeneous) is one scale factor per set, "
          "so its pooled |r| cannot exceed the ceiling row (%.2f paper subset, %.2f all sets; with run-outs included "
          "0.75 / 0.68; even one free power law per set: 0.79)." % (o["paper_subset"][0], o["all"][0]), "",
          "**The paper's r values are not comparable.** Its Table 3 within-set r differ from Table 4 recomputed "
          "(log N on log sigma_a, failures): " + ", ".join("%s %.2f vs %.2f" % (k, t3[k], abs(v)) for k, v in out["within_r"].items()) +
          ". Run-outs included, linear or semi-log variants do not reproduce them either, so the paper computes r "
          "differently or on different data (ref. [16]); its 0.84 cannot be benchmarked against this ceiling.", "",
          "**Which physics.** Shear-based Findley on the applied stresses correlates best in-sample (k = 0.1-0.2); "
          "out of sample (k re-selected per fold) its gain over the nominal stress is small (see LOSO). Normal-stress "
          "parameters do not: isotropic Walker r over all sets for gamma 0 / 0.25 / 0.5 / 0.75 / 1 = " + " / ".join(
              "%.2f" % abs(v) for _, v in out["walker_iso_scan_all"]) + ". Consistent with shear-driven matrix-fibre "
          "debonding (paper Sect. 3).", "",
          "## Leave-one-set-out prediction: rms log10 N (bias) of the held-out set; fitted k re-selected inside each fold", "",
          "| set | phase | R | " + " | ".join(LAB[n] for n in res) + " |", "|---|---|---|" + "---|" * len(res)]
    for s, m in out["sets"].items():
        L.append("| %s | %d | %g | %s |" % (s, m["delta"], m["R"], " | ".join(
            "%.2f (%+.2f)" % (res[n]["loso"][s]["rms"], res[n]["loso"][s]["bias"]) for n in res)))
    L.append("| **mean** | | | %s |" % " | ".join("**%.2f**" % np.mean([res[n]["loso"][s]["rms"] for s in out["sets"]])
                                                 for n in res))
    L += ["", "![results](cp_duda_cfrp.png)"]
    open(os.path.join(OUT, "CP_DUDA_CFRP.md"), "w").write("\n".join(L) + "\n")
    fig(out)


def fig(out):
    res = out["res"]
    N = np.array(out["N"])
    sets = np.array(out["setcol"])
    fail = np.array(out["fail"])
    colR = {-1.0: "#2e86c1", 0.05: "#e67e22", 0.5: "#8e44ad"}
    mk = {0: "o", 45: "^", 90: "s"}
    f, ax = plt.subplots(1, 5, figsize=(22, 4.4), layout="constrained")
    for a, name in zip(ax[:4], ("nominal", "paper", "cp_fibre", "cp_iso")):
        P = np.array(res[name]["P"])
        for s, m in out["sets"].items():
            sel = sets == s
            c = colR[m["R"]]
            a.loglog(N[sel & fail], P[sel & fail], mk[m["delta"]], color=c, ms=5, mfc=c, alpha=0.85,
                     label=f"{s}: {m['delta']}°, R={m['R']:g}")
            a.loglog(N[sel & ~fail], P[sel & ~fail], mk[m["delta"]], color=c, ms=5, mfc="none")
        q = res[name]["all"]
        if abs(q["r"]) >= 0.3:                                             # inverse of the N-on-P fit; meaningless if r ~ 0
            Ng = np.geomspace(N[fail].min(), N[fail].max(), 20)
            a.loglog(Ng, 10 ** ((np.log10(Ng) - q["a"]) / q["b"]), "k-", lw=1)
        a.set_ylim(P.min() / 1.15, P.max() * 1.15)
        a.set(xlabel="cycles to failure N", ylabel="parameter [MPa]",
              title=f"{LAB[name]}\nall sets: r = {q['r']:.2f}, rms log N = {q['rms']:.2f}")
        a.title.set_fontsize(9)
        a.grid(alpha=0.3, which="both")
    ax[0].legend(fontsize=6, ncol=1, loc="lower left")
    a = ax[4]
    names = list(res)
    x = np.arange(len(out["sets"]))
    w = 0.8 / len(names)
    for i, n in enumerate(names):
        a.bar(x + i * w - 0.4 + w / 2, [res[n]["loso"][s]["rms"] for s in out["sets"]], w, label=LAB[n])
    a.set_xticks(x, [f"{s}\n{m['delta']}°/{m['R']:g}" for s, m in out["sets"].items()], fontsize=7)
    a.set(ylabel="held-out rms log10 N", title="leave-one-set-out prediction")
    a.legend(fontsize=6.5)
    a.grid(alpha=0.3, axis="y")
    _save(f, OUT, "cp_duda_cfrp")


if __name__ == "__main__":
    report(main())
