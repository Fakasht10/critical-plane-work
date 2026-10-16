"""Step 1: independent Python reimplementation of the Effective Critical Plane (ECP) method of Chiocca & Frendo,
Int. J. Fatigue 189 (2024) 108565, on their own data (En3B steel, 21 notched series; see ecp_data.py).

Method (paper Sec. 2, Eqs. 5-8): FS(r_c) from stress/strain tensors averaged over a sphere (3D) or disc (2D) of
radius r_c around the hot-spot node; power law FS = c N^m fitted (log FS on log N) to two calibration series
(V-notch and Hole 8 mm, R = -1); r_c = argmin of the residual std; then N_exp = (FS / c)^(1/m) for all other series.

What this script checks:
  A. reproduction of the repository's MATLAB outputs (calibration / prevision / 19 per-series T values);
  B. the paper's "prevision T = 1.77" refits a new power law to the prevision data - here it is separated from the
     predictive error E_N = log10(N_f / N_exp) against the calibration curve (paper Eq. 9);
  C. r_c sensitivity (r_c = 0 is the hot-spot / no-averaging limit), regression direction, the MPHP corrector;
  D. robustness of r_c to the choice of calibration series (paper Table 3, 26 combinations).

    python ecp_chiocca.py          -> results/ECP_CHIOCCA.md / .json / figures
"""
import json
import os

import numpy as np

import ecp_data as ed
from plots import _save, plt

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "results")
RC = np.linspace(0.0, 0.5, 500)                                # the repository grid (MAIN.m)
D10 = 1.282                                                    # P10 / P90 quantile of the normal
GEO = {"U": "U_notch", "V": "V_notch", "H8": "Hole_8", "H3.5": "Hole_3_5"}
# paper Table 3: calibration combinations -> (r_c, slope 1/|m|)
TABLE3 = [
    ("U-1 V-1", .19, 3.62), ("U.1 V.1", .18, 3.12), ("H8-1 V-1", .20, 4.29), ("H8.1 V.1", .23, 3.39),
    ("H3.5.1 U.1", .19, 3.85), ("H3.5-1 U-1", .19, 4.39), ("H3.5-1 V-1", .19, 4.48), ("H3.5.1 V.1", .19, 3.22),
    ("H3.5-1 V-1 H3.5.1 V.1", .20, 3.69), ("H3.5-1 V-1 H3.5.1 U.1", .21, 4.18), ("H3.5-1 V-1 H3.5.1 V.1", .20, 3.69),
    ("U-1 V-1 H8-1 H3.5-1", .20, 4.10), ("U.1 V.1 H8.1 H3.5.1", .20, 3.52), ("H8-1 V-1 H8.1 V.1", .21, 3.73),
    ("H3.5-1 V-1 H8.1 V.1", .21, 3.80), ("H3.5-1 V-1 H8.1 U.1", .23, 4.42), ("U-1 V-1 U.1 V.1", .18, 3.40),
    ("U-1 V-1 U.1 H8.1", .24, 3.68), ("U-1 V-1 U.1 H3.5.1", .21, 3.62), ("H8-1 V-1 U.1 V.1", .19, 3.45),
    ("H8-1 V-1 U.1 H8.1", .23, 4.17), ("H8-1 V-1 U.1 H3.5.1", .22, 3.95), ("H3.5-1 V-1 U.1 V.1", .19, 3.57),
    ("U-1 V-1 V.1 H8.1", .20, 3.50), ("U-1 V-1 V.1 H3.5.1", .19, 3.39), ("H8-1 V-1 V.1 H3.5.1", .20, 3.57)]


# ------------------------------------------------------------------------------------------------------ statistics
def fit_cp_on_n(CP, Nf):
    """MATLAB pinv fit of log CP = log c + m log N (CP is the dependent variable), columnwise over r_c."""
    A = np.column_stack([np.ones(len(Nf)), np.log(Nf)])
    sol = np.linalg.lstsq(A, np.log(CP), rcond=None)[0]
    return np.exp(sol[0]), sol[1]


def scatter_T(CP, Nf, c, m):
    """T = c10 / c90 = exp(2 * 1.282 * std(ln N residual) * |m|) (repository FatigueScatter.m; std with ddof = 1)."""
    r = np.log(Nf)[:, None] - np.log((CP / c) ** (1 / m)) if CP.ndim == 2 else np.log(Nf) - np.log((CP / c) ** (1 / m))
    s = r.std(0, ddof=1)
    return np.exp(2 * D10 * s * np.abs(m)), s


def self_T(CP, Nf):
    """Fatigue scatter of a data set against its OWN refitted power law (what the paper's Table 4 reports)."""
    c, m = fit_cp_on_n(CP, Nf)
    return scatter_T(CP, Nf, c, m)[0], c, m


def e_n(CP, Nf, c, m):
    """Paper Eq. 9: E_N = log10(N_f / N_exp); > 0 conservative."""
    return np.log10(Nf) - np.log10((CP / c) ** (1 / m))


def summary(E):
    a = np.abs(E)
    return dict(n=int(E.size), mean=float(E.mean()), std=float(E.std(ddof=1)), rms=float(np.sqrt((E ** 2).mean())),
                in2=float((a <= np.log10(2)).mean()), in3=float((a <= np.log10(3)).mean()),
                in5=float((a <= np.log10(5)).mean()))


# ------------------------------------------------------------------------------------------------------- analyses
def stack(series, rc):
    return np.vstack([s["fs"](rc) for s in series]), np.concatenate([s["Nf"] for s in series])


def calibrate(cal, rc=RC):
    CP, Nf = stack(cal, rc)
    c, m = fit_cp_on_n(CP, Nf)
    T, s = scatter_T(CP, Nf, c, m)
    i = int(np.argmin(s))
    return dict(c=c, m=m, T=T, s=s, i_opt=i, rc_opt=float(rc[i]))


def table3(S, rc=np.linspace(0, 0.7, 701), nmin=3e4):
    """Paper Table 3 re-run: r_c of minimum scatter and 1/|m| for each calibration combination (N_f > 3e4)."""
    by = {(s["geom"], s["R"]): s for s in S if s["kind"] == "Axial"}
    rows = []
    for combo, rc_p, k_p in TABLE3:
        sel = []
        for tok in combo.split():
            g, R = (tok[:-2], -1.0) if tok.endswith("-1") else (tok[:-2], 0.1)
            s = by[(GEO[g], R)]
            keep = s["Nf"] > nmin
            sel.append(dict(s, Nf=s["Nf"][keep], fs=(lambda f, k: (lambda r: f(r)[k]))(s["fs"], keep)))
        r = calibrate(sel, rc)
        rows.append(dict(combo=combo, rc_paper=rc_p, k_paper=k_p, rc=r["rc_opt"], k=float(-1 / r["m"][r["i_opt"]])))
    return rows


def run():
    S = ed.load()
    Snc = {s["name"]: s for s in ed.load(corrector=False)}
    ref = ed.reference_results()
    cal = [s for s in S if s["role"] == "calibration"]
    prv = [s for s in S if s["role"] == "prevision"]
    out = {}

    # ---------------- A. reproduction at r_c = 0.2 (MAIN.m: val = 0.2, nearest grid point)
    C = calibrate(cal)
    ix = int(np.argmin(np.abs(RC - 0.2)))
    c, m = C["c"][ix], C["m"][ix]
    CPp, Nfp = stack(prv, RC)
    tot = self_T(CPp[:, ix], Nfp)
    rep = [("calibration T", C["T"][ix], ref["calib"][0]), ("calibration c", c, ref["calib"][2]),
           ("calibration m", m, ref["calib"][4]), ("prevision T (self-fit)", tot[0], ref["total"][0]),
           ("prevision m (self-fit)", tot[2], ref["total"][4])]
    per = []
    for s in prv:
        t = self_T(s["fs"](RC)[:, ix], s["Nf"])[0]
        per.append((s["name"], t, ref["per"][s["name"]][0]))
    resid_dev = float(np.max(np.abs(C["s"] - ref["residuals"][:, 1])))
    out["reproduction"] = dict(scalars=rep, per_series=per, residual_curve_max_abs_dev=resid_dev,
                               max_rel_dev=float(max(abs(a / b - 1) for _, a, b in rep + per)))
    out["rc_opt_repo_grid"] = C["rc_opt"]

    # ---------------- B. predictive error vs the calibration curve
    def predict(series, c, m, i):
        E = np.concatenate([e_n(s["fs"](RC)[:, i], s["Nf"], c, m) for s in series])
        return E

    groups = {"prevision (19 series)": prv,
              "flat axial, prevision": [s for s in prv if s["kind"] == "Axial"],
              "flat axial, all (paper Fig. 11a)": [s for s in S if s["kind"] == "Axial"],
              "cylindrical multiaxial + torsion (paper Fig. 11b)": [s for s in S if s["kind"] != "Axial"],
              "  in-phase (0 deg)": [s for s in S if s["phase"] == 0],
              "  out-of-phase (90 deg)": [s for s in S if s["phase"] == 90],
              "  torsion": [s for s in S if s["kind"] == "Torsion"]}
    out["predictive"] = {g: summary(predict(v, c, m, ix)) for g, v in groups.items()}
    T_cal_band = float(np.exp(2 * D10 * C["s"][ix] * abs(m)))
    Ep = predict(prv, c, m, ix)
    # predictive scatter on the paper's T scale: same formula, residuals in ln N about the CALIBRATION curve (bias kept)
    T_pred = float(np.exp(2 * D10 * np.sqrt(np.mean((Ep * np.log(10)) ** 2)) * abs(m)))
    out["T_scale"] = dict(calibration=T_cal_band, prevision_selffit=float(tot[0]),
                          prevision_predictive_rms=T_pred,
                          prevision_predictive_std=float(np.exp(2 * D10 * Ep.std(ddof=1) * np.log(10) * abs(m))))
    out["per_series_bias"] = [dict(name=s["name"], label=s["label"], n=len(s["Nf"]),
                                   **{k: v for k, v in summary(e_n(s["fs"](RC)[:, ix], s["Nf"], c, m)).items()
                                      if k in ("mean", "std")}) for s in prv]

    # ---------------- C. sensitivities
    Eall = np.array([predict(prv, C["c"][i], C["m"][i], i) for i in range(RC.size)])        # (rc, n)
    rms_rc = np.sqrt((Eall ** 2).mean(1))
    i_or = int(np.argmin(rms_rc))
    out["rc_sensitivity"] = dict(
        rc0_hotspot=summary(Eall[0]), rc02=summary(Eall[ix]), rc_opt_repo=summary(Eall[C["i_opt"]]),
        oracle=dict(rc=float(RC[i_or]), **summary(Eall[i_or])))
    # regression direction: log N on log FS (ASTM E739 style, N the dependent variable)
    CPc, Nfc = stack(cal, RC)
    b, a = np.polyfit(np.log(CPc[:, ix]), np.log(Nfc), 1)
    Er = np.log10(Nfp) - (a + b * np.log(CPp[:, ix])) / np.log(10)
    out["regression_N_on_FS"] = dict(m_equiv=float(1 / b), m_paper=float(m), **summary(Er))
    # MPHP corrector off (plain MCC shear amplitude) for the 90 deg series
    p90 = [s for s in prv if s["phase"] == 90]
    out["corrector"] = dict(
        with_MPHP=summary(predict(p90, c, m, ix)),
        MCC_only=summary(np.concatenate([e_n(Snc[s["name"]]["fs"](RC)[:, ix], s["Nf"], c, m) for s in p90])),
        factor_range=[float(min(s["corr"].min() for s in p90)), float(max(s["corr"].max() for s in p90))])

    # ---------------- D. Table 3
    t3 = table3(S)
    out["table3"] = t3
    d_rc = np.array([r["rc"] - r["rc_paper"] for r in t3])
    d_k = np.array([r["k"] - r["k_paper"] for r in t3])
    out["table3_summary"] = dict(rc_mean=float(np.mean([r["rc"] for r in t3])), rc_std=float(np.std([r["rc"] for r in t3], ddof=1)),
                                 max_abs_drc=float(np.abs(d_rc).max()), max_abs_dk=float(np.abs(d_k).max()),
                                 n_rc_match_0p01=int((np.abs(d_rc) <= 0.0105).sum()))

    figures(S, cal, prv, C, ix, Eall, rms_rc, out)
    with open(os.path.join(OUT, "ECP_CHIOCCA.json"), "w") as f:
        json.dump(out, f, indent=1, default=float)
    report(out)
    return out


# ------------------------------------------------------------------------------------------------------- figures
FAM = [("flat, axial, R = -1", "#2a78d6", lambda s: s["kind"] == "Axial" and s["R"] == -1),
       ("flat, axial, R = 0.1", "#eb6834", lambda s: s["kind"] == "Axial" and s["R"] == 0.1),
       ("cyl., torsion", "#1baf7a", lambda s: s["kind"] == "Torsion"),
       ("cyl., in-phase", "#eda100", lambda s: s["phase"] == 0),
       ("cyl., 90 deg out-of-phase", "#6250d6", lambda s: s["phase"] == 90)]


def fam(s):
    for i, (lab, col, f) in enumerate(FAM):
        if f(s):
            return i, lab, col


def figures(S, cal, prv, C, ix, Eall, rms_rc, out):
    # 1. scatter vs r_c: calibration T (paper Fig. 6c) and predictive rms of E_N on the 19 prevision series
    fig, ax = plt.subplots(1, 2, figsize=(7.2, 2.8))
    ax[0].plot(RC, C["T"], color="#2a78d6")
    ax[0].axvline(C["rc_opt"], color="#52514e", lw=1, ls="--")
    ax[0].set(xlabel="control radius $r_c$ [mm]", ylabel="calibration scatter $T_{ECP}$", ylim=(1, 2.5), xlim=(0.1, 0.5),
              title="Calibration (V-notch + hole 8, R = -1)")
    ax[0].annotate(f"min at {C['rc_opt']:.3f} mm\nT = {C['T'][C['i_opt']]:.3f}", (C["rc_opt"], C["T"][C["i_opt"]]),
                   (0.3, 2.0), arrowprops=dict(arrowstyle="-", color="#52514e", lw=0.8), fontsize=8)
    ax[1].plot(RC, rms_rc, color="#eb6834", label="rms $E_N$ (prediction)")
    ax[1].plot(RC, np.abs(Eall.mean(1)), color="#52514e", lw=1.2, ls=":", label="|mean $E_N$| (bias)")
    ax[1].axvline(C["rc_opt"], color="#52514e", lw=1, ls="--")
    ax[1].set(xlabel="control radius $r_c$ [mm]", ylabel="$E_N$ = log$_{10}$($N_f$/$N_{exp}$)", xlim=(0, 0.5),
              ylim=(0, 1.5), title="Prediction, 19 prevision series")
    ax[1].legend(fontsize=8)
    _save(fig, OUT, "ecp_scatter_vs_rc")

    # 2. FS vs N (paper Fig. 8a) with the calibration band
    c, m = C["c"][ix], C["m"][ix]
    fig, ax = plt.subplots(figsize=(5.6, 3.8))
    N = np.logspace(3.5, 6.6, 50)
    sb = C["s"][ix]
    ax.plot(N, c * N ** m, color="#0b0b0b", lw=1.5, label="calibrated $P_{50}$")
    for sg in (1, -1):
        ax.plot(N, c * (N * np.exp(sg * D10 * sb)) ** m, color="#0b0b0b", lw=0.8, ls="--",
                label="$P_{10}$-$P_{90}$" if sg == 1 else None)
    seen = set()
    for s in S:
        i, lab, col = fam(s)
        cal_ = s["role"] == "calibration"
        mk = "s" if cal_ else ("o" if s["R"] == -1 else "^")
        fc = col if (s["R"] == -1 or cal_) else "none"
        ax.scatter(s["Nf"], s["fs"](RC)[:, ix], s=22, marker=mk, facecolors=fc, edgecolors=col, lw=1.0,
                   label=(lab if lab not in seen else None), zorder=3)
        seen.add(lab)
    ax.scatter([], [], marker="s", color="#52514e", label="calibration points")
    ax.set(xscale="log", yscale="log", xlabel="cycles to failure $N_f$", ylabel=r"ECP Fatemi-Socie factor $\widetilde{FS}$",
           title=f"ECP-FS, $r_c$ = {RC[ix]:.3f} mm  (open symbols: R = 0 / 0.1)")
    ax.legend(fontsize=7, loc="lower left", ncol=2)
    _save(fig, OUT, "ecp_fs_vs_n")

    # 3. expected vs experimental life: hot spot (r_c = 0) vs ECP (r_c = 0.2)
    fig, axs = plt.subplots(1, 2, figsize=(7.4, 3.6), sharey=True)
    for ax, i, ttl in ((axs[0], 0, "hot spot, $r_c$ = 0"), (axs[1], ix, f"ECP, $r_c$ = {RC[ix]:.2f} mm")):
        ci, mi = C["c"][i], C["m"][i]
        L = np.array([1e3, 1e8])
        ax.plot(L, L, color="#0b0b0b", lw=1)
        for f_, ls in ((3, "--"), (5, ":")):
            ax.plot(L, L * f_, color="#52514e", lw=0.8, ls=ls)
            ax.plot(L, L / f_, color="#52514e", lw=0.8, ls=ls, label=f"factor {f_}" if ax is axs[0] else None)
        seen = set()
        for s in S:
            k, lab, col = fam(s)
            Ne = (s["fs"](RC)[:, i] / ci) ** (1 / mi)
            ax.scatter(s["Nf"], Ne, s=16, color=col, edgecolors="white", lw=0.4, zorder=3,
                       label=lab if (lab not in seen and ax is axs[1]) else None)
            seen.add(lab)
        E = np.concatenate([e_n(s["fs"](RC)[:, i], s["Nf"], ci, mi) for s in prv])
        sm = summary(E)
        ax.set(xscale="log", yscale="log", xlim=(1e3, 1e8), ylim=(1e3, 1e8), xlabel="experimental $N_f$",
               title=f"{ttl}\nprevision: mean $E_N$ {sm['mean']:+.2f}, std {sm['std']:.2f}, {100 * sm['in3']:.0f} % in x3")
        ax.set_aspect("equal")
    axs[0].set_ylabel("expected $N_{exp}$")
    axs[0].legend(fontsize=7, loc="upper left")
    axs[1].legend(fontsize=6.5, loc="lower right")
    _save(fig, OUT, "ecp_life_hotspot_vs_ecp")

    # 4. per-series bias (the information a pooled scatter index hides)
    B = out["per_series_bias"]
    order = sorted(range(len(B)), key=lambda j: (fam(prv[j])[0], B[j]["mean"]))
    fig, ax = plt.subplots(figsize=(6.0, 4.6))
    for y, j in enumerate(order):
        col = fam(prv[j])[2]
        ax.errorbar(B[j]["mean"], y, xerr=B[j]["std"], fmt="o", color=col, ms=5, lw=1.2, capsize=0)
    ax.axvline(0, color="#0b0b0b", lw=1)
    for f_ in (3,):
        for sg in (1, -1):
            ax.axvline(sg * np.log10(f_), color="#52514e", lw=0.8, ls="--")
    ax.set_yticks(range(len(order)), [B[j]["label"] for j in order], fontsize=7)
    ax.set(xlabel="$E_N$ = log$_{10}$($N_f$ / $N_{exp}$)   (mean ± std per series; > 0 conservative)",
           title=f"Per-series prediction bias, ECP-FS $r_c$ = {RC[ix]:.2f} mm (dashed: factor 3)")
    ax.grid(axis="y", visible=False)
    _save(fig, OUT, "ecp_series_bias")


# ------------------------------------------------------------------------------------------------------- report
def report(o):
    R, P, T = o["reproduction"], o["predictive"], o["T_scale"]
    f = lambda d: f"{d['mean']:+.3f} | {d['std']:.3f} | {100 * d['in2']:.0f} | {100 * d['in3']:.0f} | {100 * d['in5']:.0f}"
    L = ["# ECP (Chiocca & Frendo 2024) reimplemented and stress-tested on its own data",
         "",
         "**Summary.**",
         f"* **Reproduction is exact.** The Python port reproduces the repository's MATLAB outputs to "
         f"{R['max_rel_dev']:.0e} relative: calibration T = {R['scalars'][0][1]:.4f}, prevision T = "
         f"{R['scalars'][3][1]:.4f} and all 19 per-series T values. The calibration residual curve matches to "
         f"{R['residual_curve_max_abs_dev']:.0e}.",
         f"* **The headline prevision T = {T['prevision_selffit']:.2f} is a self-fit, but here it happens to be "
         f"honest.** Like the paper's Table 4, it refits a new power law (m = {R['scalars'][4][1]:.3f}, against the "
         f"calibrated m = {R['scalars'][2][1]:.3f}). Scored against the calibrated curve, as a prediction must be, "
         f"the same T formula gives {T['prevision_predictive_rms']:.2f}.",
         f"* **The near-zero multiaxial bias is a cancellation.** In-phase series are non-conservative (mean E_N "
         f"{P['  in-phase (0 deg)']['mean']:+.2f}) and 90 deg series conservative "
         f"({P['  out-of-phase (90 deg)']['mean']:+.2f}). The 90 deg result rests on the MPHP corrector: with the "
         f"plain MCC shear amplitude (what `cplane` and most codes use) the 90 deg mean is "
         f"{o['corrector']['MCC_only']['mean']:+.2f}, about a factor {10 ** -o['corrector']['MCC_only']['mean']:.0f} "
         f"non-conservative.",
         f"* **Predicted lives are still good for a one-parameter model.** On the 19 prevision series, E_N has mean "
         f"{P['prevision (19 series)']['mean']:+.2f} and std {P['prevision (19 series)']['std']:.2f}; "
         f"{100 * P['prevision (19 series)']['in3']:.0f} % fall within a factor of 3 "
         f"(the paper reports 74 % on all data).",
         f"* **Averaging is what makes it work.** At r_c = 0 (hot spot, no averaging) the rms E_N is "
         f"{o['rc_sensitivity']['rc0_hotspot']['rms']:.2f} decades, against "
         f"{o['rc_sensitivity']['rc02']['rms']:.2f} at r_c = 0.2 mm. The oracle r_c, chosen on the prevision data, "
         f"is {o['rc_sensitivity']['oracle']['rc']:.3f} mm with rms {o['rc_sensitivity']['oracle']['rms']:.2f}. So "
         f"the two-series calibration picks a near-optimal radius.",
         f"* **r_c is robust to the calibration choice.** Re-running the paper's 26 combinations (Table 3) gives "
         f"r_c = {o['table3_summary']['rc_mean']:.3f} ± {o['table3_summary']['rc_std']:.3f} mm "
         f"(paper: 0.204 ± 0.016); {o['table3_summary']['n_rc_match_0p01']} of 26 radii match the paper within 0.01 mm.",
         "",
         "## A. Reproduction of the repository outputs (r_c = 0.2 mm grid point)",
         "| quantity | Python | repository (MATLAB) |", "|---|---|---|"]
    L += [f"| {n} | {a:.6g} | {b:.6g} |" for n, a, b in R["scalars"]]
    L += ["", "| prevision series | T (self-fit), Python | repository |", "|---|---|---|"]
    L += [f"| {n} | {a:.5f} | {b:.5f} |" for n, a, b in R["per_series"]]
    L += [f"", f"Optimal r_c on the repository grid (0-0.5 mm, 500 points): {o['rc_opt_repo_grid']:.4f} mm. MAIN.m "
          f"evaluates at the fixed value 0.2.", "",
          "## B. Predictive error against the calibrated curve, r_c = 0.2 mm",
          "E_N = log10(N_f / N_exp), paper Eq. 9; > 0 means conservative. All rows are predictions except where a "
          "group includes the two calibration series (marked 'all').", "",
          "| group | n | mean E_N | std E_N | % in x2 | % in x3 | % in x5 |", "|---|---|---|---|---|---|---|"]
    L += [f"| {g} | {d['n']} | {f(d)} |" for g, d in P.items()]
    L += ["", "Paper Fig. 11 (FS ECP): uniaxial mean -0.022, std 0.285; multiaxial mean +0.013, std 0.492.",
          "",
          "**Scatter on the paper's T scale** (T = exp(2 · 1.282 · s · |m|), s = std of ln N):", "",
          "| definition | T |", "|---|---|",
          f"| calibration data vs calibration curve (paper: 1.20) | {T['calibration']:.3f} |",
          f"| prevision vs its own refitted curve (paper Table 4: 1.77) | {T['prevision_selffit']:.3f} |",
          f"| prevision vs calibration curve, bias removed | {T['prevision_predictive_std']:.3f} |",
          f"| prevision vs calibration curve, rms (bias kept) | {T['prevision_predictive_rms']:.3f} |",
          "", "Per-series bias (mean ± std of E_N) is in `ecp_series_bias.png`. Every series has a T near 1.1-1.4 "
          "about its own line (Table 4), but the series lines are offset from each other. The offsets are what a "
          "designer would actually incur.", "",
          "## C. Sensitivities (19 prevision series)",
          "| variant | mean E_N | std | % x2 | % x3 | % x5 |", "|---|---|---|---|---|---|"]
    rs = o["rc_sensitivity"]
    L += [f"| r_c = 0 (hot spot, no averaging) | {f(rs['rc0_hotspot'])} |",
          f"| r_c = 0.2 mm (paper) | {f(rs['rc02'])} |",
          f"| r_c = {o['rc_opt_repo_grid']:.3f} mm (calibration optimum) | {f(rs['rc_opt_repo'])} |",
          f"| r_c = {rs['oracle']['rc']:.3f} mm (oracle on prevision data; not a prediction) | {f(rs['oracle'])} |",
          f"| r_c = 0.2, regression ln N on ln FS (equivalent m = {o['regression_N_on_FS']['m_equiv']:.3f} vs "
          f"{o['regression_N_on_FS']['m_paper']:.3f}) | {f(o['regression_N_on_FS'])} |",
          "", "90 deg out-of-phase series only (the corrector multiplies FS by the MPHP / MCC ratio, "
          f"{o['corrector']['factor_range'][0]:.2f}-{o['corrector']['factor_range'][1]:.2f}):", "",
          "| shear amplitude | mean E_N | std | % x2 | % x3 | % x5 |", "|---|---|---|---|---|---|",
          f"| MPHP (paper) | {f(o['corrector']['with_MPHP'])} |",
          f"| MCC (no corrector) | {f(o['corrector']['MCC_only'])} |", "",
          "## D. Paper Table 3 re-run (N_f > 3e4, r_c grid 0-0.7 mm in 0.001 steps)",
          "| combination (geometry, R) | r_c paper | r_c here | 1/|m| paper | 1/|m| here |", "|---|---|---|---|---|"]
    L += [f"| {r['combo']} | {r['rc_paper']:.2f} | {r['rc']:.3f} | {r['k_paper']:.2f} | {r['k']:.2f} |" for r in o["table3"]]
    L += ["", f"Max |Δr_c| = {o['table3_summary']['max_abs_drc']:.3f} mm, max |Δ(1/m)| = {o['table3_summary']['max_abs_dk']:.2f}. "
          "Table rows 9 and 11 are identical in the paper; both are kept as printed.", "",
          "## Theory / physics / numerics / visualisation checklist",
          "* **Theory.** Fatemi-Socie on the averaged tensors: FS = (Δγ/2)(1 + k σ_n,max / σ_y), k = 0.4, "
          "σ_y = 606.2 MPa. Linear elasticity, so FS(F) = (F Δγ₁ / 2)(1 + k F σ_n1 / σ_y) from unit-load fields.",
          "* **Physics.** En3B is a ductile low-carbon steel with shear-dominated crack initiation, which is why FS "
          "suits it (the paper's SWT variant scatters 1:4.55). Only the critical point and the averaging radius "
          "carry the size effect; plasticity at short lives is ignored (the calibration keeps N_f > 3e4).",
          "* **Numerics.** The authors' degree-6 polynomial fits of the averaged quantities vs r_c are used as given "
          "(R² > 99.9 %). Averaging is a plain node mean in a sphere or disc on their ANSYS meshes (0.01-0.1 mm). "
          "Step 2 (`ecp_fe.py`) recomputes these factors with our NGSolve models and `cplane`.",
          "* **Visualisation.** Calibration and prediction are reported separately throughout. Figures: "
          "`ecp_scatter_vs_rc`, `ecp_fs_vs_n`, `ecp_life_hotspot_vs_ecp`, `ecp_series_bias`.",
          "", "Data: github.com/achiocca1/ECP @ 9d51152 (GPL-3, fetched by `ecp_data.fetch()`, not redistributed)."]
    with open(os.path.join(OUT, "ECP_CHIOCCA.md"), "w") as fh:
        fh.write("\n".join(L) + "\n")


if __name__ == "__main__":
    o = run()
    print(json.dumps({k: o[k] for k in ("T_scale", "rc_sensitivity", "regression_N_on_FS", "corrector", "table3_summary")},
                     indent=1, default=float))
    print(json.dumps(o["predictive"], indent=1))
