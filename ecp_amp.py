"""Shear-amplitude measure for non-proportional loading: which one centres the 90 deg out-of-phase ECP series?

Data: Chiocca & Frendo (2024) / Susmel & Taylor (2008) En3B V-notched bars, tension-torsion, 90 deg out of phase,
R = -1 and 0, root radii 0.2 / 1.25 / 4 mm (6 series, 41 specimens). Factors: our NGSolve + cplane ECP factors
(ecp_fe.py, cached), r_c = 0.2 mm, calibration on V-notch + hole 8 at R = -1 (proportional, so it is identical for
every measure). Only the shear amplitude of the 90 deg paths changes:
  MCC     minimum circumscribed circle (cplane default; Papadopoulos / Dang Van school);
  MRH     maximum rectangular hull, Araujo et al. 2011;  ellipse a, b -> sqrt(a^2 + b^2);
  length  path length / 4 (a half-perimeter measure in the spirit of their MPHP);  circle -> pi a / 2;
  their MPHP corrector (factor per specimen from their repository) on our MCC factors, for reference.
No parameter is fitted to the 90 deg data: every row is a prediction.

    python ecp_amp.py   -> results/ECP_AMPLITUDE.md / .json / ecp_amplitude.png
"""
import json
import os

import numpy as np

import cplane as cp
import ecp_chiocca as ec
import ecp_data as ed
import ecp_fe as F
from plots import _save, plt

OUT = F.OUT
MEASURES = [("MCC", "mcc"), ("MRH", "mrh"), ("path length / 4", "length")]


def run():
    bars, plates = F.factors()
    S_nc = ed.load(corrector=False)
    i2 = int(np.argmin(np.abs(F.RCS - 0.2)))
    variants = {}
    for lab, amp in MEASURES:
        variants[lab] = F.our_series(S_nc, bars, plates, corrector=False, amplitude=amp)
    variants["MCC x their MPHP corrector"] = F.our_series(ed.load(), bars, plates, corrector=True)
    variants["their factors + MPHP (paper)"] = ed.load()

    out = dict(variants={}, per_series={}, ratio={})
    for tag, ser in variants.items():
        cal = [s for s in ser if s["role"] == "calibration"]
        prv = [s for s in ser if s["role"] == "prevision"]
        C = ec.calibrate(cal, F.RCS)
        c, m = C["c"][i2], C["m"][i2]
        E = lambda grp: np.concatenate([ec.e_n(s["fs"](F.RCS)[:, i2], s["Nf"], c, m) for s in grp])
        p90 = [s for s in prv if s["phase"] == 90]
        out["variants"][tag] = dict(m=float(m), c=float(c), all=ec.summary(E(prv)), outphase=ec.summary(E(p90)),
                                    inphase=ec.summary(E([s for s in prv if s["phase"] == 0])),
                                    p90_R1=ec.summary(E([s for s in p90 if s["R"] == -1])),
                                    p90_R0=ec.summary(E([s for s in p90 if s["R"] == 0])))
        out["per_series"][tag] = {s["label"]: ec.summary(E([s])) for s in p90}

    # amplitude ratio (FS_measure / FS_MCC) per 90 deg specimen vs their corrector
    mcc = {s["name"]: s for s in variants["MCC"] if s["phase"] == 90}
    for lab, _ in MEASURES[1:]:
        out["ratio"][lab] = {s["label"]: float(np.mean(s["fs_tab"][:, i2] / mcc[s["name"]]["fs_tab"][:, i2]))
                             for s in variants[lab] if s["phase"] == 90}
    out["ratio"]["their MPHP corrector"] = {s["label"]: float(s["corr"].mean()) for s in ed.load() if s["phase"] == 90}

    # shear path on the critical plane, first specimen of V-notch r = 1.25, 90 deg, R = -1 (for the figure)
    s = next(x for x in S_nc if x["name"] == "Multi_V_notch_12_Ph_90_R_-1")
    b = bars["V_notch_12"]
    t = np.linspace(0, 2 * np.pi, 360, endpoint=False)
    L = np.column_stack([s["S"][0] * np.sin(t), s["tau"][0] * np.sin(t - np.pi / 2)])
    H = cp.superpose(np.stack([b["A_ax"][i2:i2 + 1], b["A_to"][i2:i2 + 1]]), L)
    r = cp.search_path(H, F.FS)
    Y = H[0] @ cp.proj_rows(r["normal"][0]).T
    out["example"] = dict(name=s["label"], sa=float(s["S"][0]), ta=float(s["tau"][0]), normal=r["normal"][0].tolist())
    out["calibration_identical"] = bool(len({round(v["m"], 12) for k, v in out["variants"].items()
                                             if k != "their factors + MPHP (paper)"}) == 1)
    with open(os.path.join(OUT, "ECP_AMPLITUDE.json"), "w") as f:
        json.dump(out, f, indent=1, default=float)
    figure(out, Y)
    report(out)
    return out


COL = {"MCC": "#2a78d6", "MRH": "#1baf7a", "path length / 4": "#eda100", "MCC x their MPHP corrector": "#6250d6",
       "their factors + MPHP (paper)": "#52514e"}


def figure(o, Y):
    fig, ax = plt.subplots(1, 2, figsize=(9.4, 3.8), gridspec_kw=dict(width_ratios=[1, 1.7]))
    u, v = Y[:, 1] - Y[:, 1].mean(), Y[:, 2] - Y[:, 2].mean()
    ax[0].plot(np.r_[u, u[:1]], np.r_[v, v[:1]], color="#0b0b0b", lw=1.6, label="shear path")
    R = cp.shear_amplitude(u[None], v[None], "mcc")[0]
    th = np.linspace(0, 2 * np.pi, 200)
    ax[0].plot(R * np.cos(th), R * np.sin(th), color=COL["MCC"], lw=1.2, label=f"MCC: {R:.0f}")
    Rm = cp.shear_amplitude(u[None], v[None], "mrh")[0]
    ax[0].plot(Rm * np.cos(th), Rm * np.sin(th), color=COL["MRH"], lw=1.0, ls="--", label=f"MRH: {Rm:.0f}")
    Rl = cp.shear_amplitude(u[None], v[None], "length")[0]
    ax[0].plot(Rl * np.cos(th), Rl * np.sin(th), color=COL["path length / 4"], lw=1.0, ls=":",
               label=f"length / 4: {Rl:.0f}")
    ax[0].set_aspect("equal")
    ax[0].set(xlabel=r"$\tau_u$ [MPa]", ylabel=r"$\tau_v$ [MPa]",
              title=f"Critical-plane shear path\n{o['example']['name']}, " + "$r_c$ = 0.2 mm", )
    ax[0].legend(fontsize=7, loc="lower left")
    tags = list(o["per_series"])
    labels = list(o["per_series"][tags[0]])
    for k, tag in enumerate(tags):
        y = np.arange(len(labels)) + (k - (len(tags) - 1) / 2) * 0.15
        mu = [o["per_series"][tag][l]["mean"] for l in labels]
        sd = [o["per_series"][tag][l]["std"] for l in labels]
        ax[1].errorbar(mu, y, xerr=sd, fmt="o", ms=4, lw=1, color=COL[tag], label=tag)
    ax[1].axvline(0, color="#0b0b0b", lw=1)
    for sg in (1, -1):
        ax[1].axvline(sg * np.log10(3), color="#52514e", lw=0.8, ls="--")
    ax[1].set_yticks(range(len(labels)), [l.replace(", 90 deg", "") for l in labels], fontsize=8)
    ax[1].invert_yaxis()
    ax[1].grid(axis="y", visible=False)
    ax[1].set(xlabel="$E_N$ = log$_{10}$($N_f/N_{exp}$), mean ± std  (> 0 conservative; dashed: factor 3)",
              title="90 deg out-of-phase series: prediction error by shear-amplitude measure")
    ax[1].legend(fontsize=7, loc="upper left", bbox_to_anchor=(1.01, 1.0))
    fig.tight_layout()
    _save(fig, OUT, "ecp_amplitude")


def report(o):
    V = o["variants"]
    f = lambda d: f"{d['mean']:+.2f} | {d['std']:.2f} | {100 * d['in2']:.0f} | {100 * d['in3']:.0f}"
    best = min(["MCC", "MRH", "path length / 4"], key=lambda k: V[k]["outphase"]["rms"])
    L = ["# Shear-amplitude measure for 90 deg out-of-phase loading (ECP, En3B V-notched bars)", "",
         "**Summary.**",
         f"* **Best measure: {best}.** On the 41 out-of-phase specimens (6 series), the mean prediction error E_N is:",
         f"  * MCC: {V['MCC']['outphase']['mean']:+.2f} (lives {10 ** -V['MCC']['outphase']['mean']:.1f}x too "
         "long, unsafe).",
         f"  * MRH: {V['MRH']['outphase']['mean']:+.2f}.",
         f"  * path length / 4: {V['path length / 4']['outphase']['mean']:+.2f}.",
         f"  * their MPHP corrector: {V['MCC x their MPHP corrector']['outphase']['mean']:+.2f}.",
         f"  {best} is a prediction: nothing was fitted to these data. Calibration and all proportional series are "
         f"unchanged ({'identical calibration for every measure' if o['calibration_identical'] else 'CHECK: calibration differs'}).",
         f"* **Within a factor of 3:** {100 * V[best]['outphase']['in3']:.0f} % with {best}, against "
         f"{100 * V['MCC']['outphase']['in3']:.0f} % with MCC and "
         f"{100 * V['MCC x their MPHP corrector']['outphase']['in3']:.0f} % with their corrector.",
         f"* **Remaining split by R** ({best}): R = -1 mean {V[best]['p90_R1']['mean']:+.2f}, R = 0 mean "
         f"{V[best]['p90_R0']['mean']:+.2f}.",
         "* **Physics.** Out-of-phase loading of a ductile steel causes extra cyclic (non-proportional) hardening and "
         "more active slip systems. MCC sees only the size of the shear path, not how round it is. MRH and the "
         "path length grow with roundness (by up to sqrt(2) and pi/2 for a circular path), which is the direction "
         "the data require. A material-specific non-proportionality constant (Kanazawa / Itoh type) would be the "
         "next refinement, but it needs an extra calibration series.",
         "",
         "## All 19 prevision series and by family (r_c = 0.2 mm)",
         "Cells: mean E_N | std | % within factor 2 | % within factor 3.", "",
         "| variant | all 19 | in-phase | 90 deg | 90 deg, R = -1 | 90 deg, R = 0 |", "|---|---|---|---|---|---|"]
    for tag, d in V.items():
        L.append(f"| {tag} | {f(d['all'])} | {f(d['inphase'])} | {f(d['outphase'])} | {f(d['p90_R1'])} | {f(d['p90_R0'])} |")
    L += ["", "## Per 90 deg series", "| series | " + " | ".join(V) + " |", "|---|" + "---|" * len(V)]
    for lab in o["per_series"]["MCC"]:
        L.append(f"| {lab} | " + " | ".join(f"{o['per_series'][t][lab]['mean']:+.2f} ± {o['per_series'][t][lab]['std']:.2f}"
                                         for t in V) + " |")
    L += ["", "## Amplitude ratio FS(measure) / FS(MCC) at r_c = 0.2 mm, mean over specimens",
          "| series | MRH | path length / 4 | their MPHP corrector |", "|---|---|---|---|"]
    for lab in o["ratio"]["MRH"]:
        L.append(f"| {lab} | {o['ratio']['MRH'][lab]:.3f} | {o['ratio']['path length / 4'][lab]:.3f} | "
                 f"{o['ratio']['their MPHP corrector'][lab]:.3f} |")
    L += ["", "Their corrector (1.46-1.77) exceeds pi/2 = 1.571, the largest possible length/MCC ratio on a single "
          "plane, in 4 of 6 series. So their MPHP is not simply the half-perimeter on the MCC plane; the paper does "
          "not define it fully.", "",
          "## Checklist",
          "* **Theory.** FS = (τ_a / G)(1 + k σ_n,max / σ_y) on the plane that maximises it, separately for each "
          "measure. All three measures equal the half-range on a straight (proportional) path, so only the 90 deg "
          "series can change.",
          "* **Physics.** Linear elastic, with no plasticity correction. Non-proportional hardening enters only "
          "through the shape of the path.",
          "* **Numerics.**",
          "  * Path: 72 samples per cycle.",
          "  * MRH: 48-angle scan plus golden-section polish.",
          "  * Path length: the closed polygon.",
          "  * Kernel: verified against closed forms (ellipse, line) and a 30k-normal brute force "
          "(`tests/test_cplane.py`).",
          "  * Screening bounds: √2·Tresca (MRH) and the Tresca path length (length).",
          "* **Visualisation.** `ecp_amplitude.png`:",
          "  * left: the critical-plane shear path for one specimen, with the three amplitude radii;",
          "  * right: the per-series prediction error for each variant.",
          ""]
    with open(os.path.join(OUT, "ECP_AMPLITUDE.md"), "w") as fh:
        fh.write("\n".join(L) + "\n")


if __name__ == "__main__":
    o = run()
    print(json.dumps({k: {g: o["variants"][k][g] for g in ("outphase", "p90_R1", "p90_R0")} for k in o["variants"]},
                     indent=1))
