"""Volume average (ECP) vs point method (TCD PM) on the Chiocca & Frendo / Susmel & Taylor En3B data, with the
Fatemi-Socie and Smith-Watson-Topper critical-plane parameters - all from our NGSolve fields and cplane.

Same calibration for every method (paper's choice): V-notch + hole 8 mm plates at R = -1 (12 tests, N_f > 3e4);
the length parameter x (control radius r_c for ECP, bisector distance d for PM) is the one minimising the
calibration scatter.  x(N) = A (N / 1e5)^B is the life-dependent variant (Susmel & Taylor's L(N) for PM, the
repository's 'VAR' option for ECP).  Every other series (19) is a prediction.

Parameters (linear elastic, En3B: E = 200 GPa, nu = 0.3, sigma_y = 606.2 MPa, k = 0.4):
  FS  = gamma_a (1 + k sn_max / sigma_y) on the plane maximising it; gamma_a = tau_a / G with tau_a by MCC or by
        path length / 4 (ecp_amp.py);
  SWT = eps_n,a * sn_max on the plane of maximum normal-strain amplitude (paper Eq. 2, Socie's form).

    python ecp_compare.py   -> results/ECP_COMPARE.md / .json / ecp_compare_*.png
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
X = F.RCS
NT = 36
NREF = 1e5
# paper Table 4, SWT column (r_c = 0.25 mm), and the combined values
PAPER_SWT = {"Axial_U_notch_R_-1": 1.111, "Axial_V_notch_R_-1": 1.225, "Axial_Hole_3_5_R_-1": 1.101,
             "Axial_Hole_8_R_-1": 1.106, "Torsion_V_notch_02_R_-1": 1.426, "Multi_V_notch_02_Ph_0_R_-1": 1.522,
             "Multi_V_notch_02_Ph_90_R_-1": 1.650, "Multi_V_notch_12_Ph_0_R_-1": 1.178,
             "Multi_V_notch_12_Ph_90_R_-1": 1.364, "Multi_V_notch_4_Ph_0_R_-1": 1.385,
             "Multi_V_notch_4_Ph_90_R_-1": 1.330, "Axial_U_notch_R_0.1": 1.156, "Axial_V_notch_R_0.1": 1.598,
             "Axial_Hole_3_5_R_0.1": 1.266, "Axial_Hole_8_R_0.1": 1.187, "Multi_V_notch_02_Ph_0_R_0": 1.217,
             "Multi_V_notch_02_Ph_90_R_0": 1.508, "Multi_V_notch_12_Ph_0_R_0": 1.108,
             "Multi_V_notch_12_Ph_90_R_0": 1.355, "Multi_V_notch_4_Ph_0_R_0": 1.468,
             "Multi_V_notch_4_Ph_90_R_0": 1.162}
PAPER = dict(SWT=dict(rc=0.25, T_cal=1.318, T_prev=4.549, uni=(0.304, 0.264), multi=(-0.494, 0.653)),
             FS=dict(rc=0.20, T_cal=1.202, T_prev=1.773, uni=(-0.0221, 0.285), multi=(0.0128, 0.492)))


# ------------------------------------------------------------------------------------------------ histories
def histories(s, T):
    """(n_specimens, n_x, NT, 6) stress history over one cycle at every x, from unit tensors T (dict)."""
    t = np.linspace(0, 2 * np.pi, NT, endpoint=False)
    if s["kind"] == "Axial":                                     # load column = F_max (R = -1: amplitude)
        f = (1 + s["R"]) / 2 + (1 - s["R"]) / 2 * np.sin(t)
        return s["load"][:, None, None, None] * f[None, None, :, None] * T["flat"][None, :, None, :]
    if s["kind"] == "Torsion":                                   # column 2 = nominal net tau_a
        return s["S"][:, None, None, None] * np.sin(t)[None, None, :, None] * T["to"][None, :, None, :]
    mr = 0.0 if s["R"] == -1 else 1.0
    sig = s["S"][:, None] * (mr + np.sin(t))[None]
    tau = s["tau"][:, None] * (mr + np.sin(t - np.deg2rad(s["phase"])))[None]
    return sig[:, None, :, None] * T["ax"][None, :, None, :] + tau[:, None, :, None] * T["to"][None, :, None, :]


def fs_values(H, amplitude):
    n, nx = H.shape[:2]
    r = cp.search_path(H.reshape(n * nx, NT, 6), F.FS, amplitude=amplitude)
    return r["value"].reshape(n, nx) / F.G


def _cap(n0, half, k):
    u, v = cp.tangent_basis(n0)
    g = np.linspace(-half, half, k)
    a, b = np.meshgrid(g, g)
    N = n0[None] + a.ravel()[:, None] * u[None] + b.ravel()[:, None] * v[None]
    return N / np.linalg.norm(N, axis=1, keepdims=True)


N_SWT = cp.fib_hemisphere(4000)


def swt_values(H):
    """eps_n,a * sn_max on the plane of max normal-strain amplitude (4000-normal grid + two 21x21 cap refinements)."""
    n, nx = H.shape[:2]
    Hf = H.reshape(-1, NT, 6)
    tr = Hf[..., :3].sum(-1, keepdims=True)
    Ef = ((1 + F.NU) * Hf - F.NU * tr * np.r_[1, 1, 1, 0, 0, 0]) / F.E      # tensor (not engineering) shear strains
    out = np.empty(len(Hf))
    for i in range(len(Hf)):
        Nn = N_SWT
        for half, k in ((None, None), (0.06, 21), (0.006, 21)):
            if half is not None:
                Nn = _cap(best, half, k)
            row = cp.proj_rows(Nn)[:, 0, :]
            en = Ef[i] @ row.T                                         # (NT, M)
            ea = 0.5 * (en.max(0) - en.min(0))
            j = int(np.argmax(ea))
            best = Nn[j]
        sn = Hf[i] @ cp.proj_rows(best[None])[0, 0]
        out[i] = ea[j] * sn.max()
    return out.reshape(n, nx)


# ------------------------------------------------------------------------------------------------ calibration
def series_with(S, values):
    out = []
    for s in S:
        s = dict(s)
        tab = values[s["name"]]
        s["tab"] = tab
        s["fs"] = (lambda tb: (lambda x: np.array([np.interp(x, X, r) for r in tb])))(tab)
        out.append(s)
    return out


def const_method(ser):
    cal = [s for s in ser if s["role"] == "calibration"]
    prv = [s for s in ser if s["role"] == "prevision"]
    C = ec.calibrate(cal, X)
    i = C["i_opt"]
    c, m = C["c"][i], C["m"][i]
    E = {s["name"]: ec.e_n(s["tab"][:, i], s["Nf"], c, m) for s in prv}
    return dict(x=float(X[i]), T_cal=float(C["T"][i]), c=float(c), m=float(m), E=E, T_curve=C["T"].tolist(),
                T_prev_self=float(ec.self_T(np.concatenate([s["tab"][:, i] for s in prv]),
                                            np.concatenate([s["Nf"] for s in prv]))[0]),
                T_series_self={s["name"]: float(ec.self_T(s["tab"][:, i], s["Nf"])[0]) for s in ser})


def _at(tab, x):
    """tab (n, len(X)) at per-specimen x (n,)."""
    return np.array([np.interp(xi, X, r) for xi, r in zip(x, tab)])


def var_method(ser, A_grid=np.linspace(0.01, 0.5, 50), B_grid=np.linspace(-0.6, 0.3, 37)):
    cal = [s for s in ser if s["role"] == "calibration"]
    prv = [s for s in ser if s["role"] == "prevision"]
    Nf = np.concatenate([s["Nf"] for s in cal])
    tabs = np.vstack([s["tab"] for s in cal])
    best = (np.inf, None)
    for A in A_grid:
        for B in B_grid:
            x = np.clip(A * (Nf / NREF) ** B, 0, X[-1])
            P = _at(tabs, x)
            if np.any(P <= 0):
                continue
            c, m = ec.fit_cp_on_n(P[:, None], Nf)
            s_ = np.std(np.log(Nf) - np.log((P / c[0]) ** (1 / m[0])), ddof=1)
            T = np.exp(2 * ec.D10 * s_ * abs(m[0]))
            if T < best[0]:
                best = (T, (A, B, c[0], m[0]))
    T_cal, (A, B, c, m) = best
    Ng = np.logspace(2.5, 8, 600)
    xg = np.clip(A * (Ng / NREF) ** B, 0, X[-1])
    E = {}
    for s in prv:
        e = []
        for row, nf in zip(s["tab"], s["Nf"]):
            f = np.log(np.interp(xg, X, row)) - np.log(c) - m * np.log(Ng)
            k = np.where(np.diff(np.sign(f)) != 0)[0]
            if len(k):                                                  # first crossing (shortest life)
                k = k[0]
                Ne = np.exp(np.interp(0, [f[k + 1], f[k]][::1 if f[k + 1] < f[k] else -1],
                                      [np.log(Ng[k + 1]), np.log(Ng[k])][::1 if f[k + 1] < f[k] else -1]))
            else:
                Ne = Ng[-1] if f[-1] < 0 else Ng[0]
            e.append(np.log10(nf / Ne))
        E[s["name"]] = np.array(e)
    return dict(A=float(A), B=float(B), T_cal=float(T_cal), c=float(c), m=float(m), E=E,
                x_at=dict(N1e4=float(np.clip(A * 0.1 ** B, 0, X[-1])), N1e6=float(np.clip(A * 10 ** B, 0, X[-1]))))


GROUPS = {"all 19": lambda s: True, "flat axial": lambda s: s["kind"] == "Axial",
          "torsion": lambda s: s["kind"] == "Torsion", "in-phase": lambda s: s["phase"] == 0,
          "90 deg": lambda s: s["phase"] == 90}


def grouped(res, S):
    prv = {s["name"]: s for s in S if s["role"] == "prevision"}
    return {g: ec.summary(np.concatenate([e for n, e in res["E"].items() if f(prv[n])])) for g, f in GROUPS.items()}


# ------------------------------------------------------------------------------------------------ run
def run():
    bars, plates = F.factors()
    S = ed.load(corrector=False)
    fields = {}
    for kind, key in (("ECP", "A"), ("PM", "L")):
        fields[kind] = {}
        for s in S:
            g = s["geom"]
            if s["kind"] == "Axial":
                T = dict(flat=plates[g][key])
            else:
                b = bars[g]
                T = dict(ax=b[key + "_ax"] if key == "L" else b["A_ax"], to=b[key + "_to"] if key == "L" else b["A_to"])
            fields[kind][s["name"]] = histories(s, T)
    params = {"FS (MCC)": lambda H: fs_values(H, "mcc"), "FS (length/4)": lambda H: fs_values(H, "length"),
              "SWT": swt_values}
    out = dict(methods={}, paper=PAPER)
    for pname, fn in params.items():
        for kind in ("ECP", "PM"):
            vals = {n: fn(H) for n, H in fields[kind].items()}
            ser = series_with(S, vals)
            cst, var = const_method(ser), var_method(ser)
            tag = f"{pname}, {kind}"
            out["methods"][tag] = dict(const=dict({k: v for k, v in cst.items() if k != "E"}, groups=grouped(cst, S),
                                                  per={n: float(e.mean()) for n, e in cst["E"].items()}),
                                       var=dict({k: v for k, v in var.items() if k != "E"}, groups=grouped(var, S)))
            if pname == "SWT" and kind == "ECP":                       # paper's own r_c = 0.25 for SWT
                i = int(np.argmin(np.abs(X - 0.25)))
                cal = [s for s in ser if s["role"] == "calibration"]
                C = ec.calibrate(cal, X)
                c, m = C["c"][i], C["m"][i]
                uni = np.concatenate([ec.e_n(s["tab"][:, i], s["Nf"], c, m) for s in ser if s["kind"] == "Axial"])
                mul = np.concatenate([ec.e_n(s["tab"][:, i], s["Nf"], c, m) for s in ser if s["kind"] != "Axial"])
                prv = [s for s in ser if s["role"] == "prevision"]
                out["swt_repro"] = dict(
                    T_cal=float(C["T"][i]),
                    T_prev=float(ec.self_T(np.concatenate([s["tab"][:, i] for s in prv]),
                                           np.concatenate([s["Nf"] for s in prv]))[0]),
                    uni=ec.summary(uni), multi=ec.summary(mul),
                    per={s["name"]: float(ec.self_T(s["tab"][:, i], s["Nf"])[0]) for s in ser})
            print(tag, f"x = {cst['x']:.3f}", out["methods"][tag]["const"]["groups"]["all 19"], flush=True)
    with open(os.path.join(OUT, "ECP_COMPARE.json"), "w") as f:
        json.dump(out, f, indent=1, default=float)
    figures(out)
    report(out)
    return out


COL = {"FS (MCC)": "#2a78d6", "FS (length/4)": "#1baf7a", "SWT": "#eb6834"}


def figures(o):
    fig, ax = plt.subplots(1, 2, figsize=(8.6, 3.0), sharey=True)
    for k, kind in enumerate(("ECP", "PM")):
        for p in COL:
            m = o["methods"][f"{p}, {kind}"]["const"]
            ax[k].plot(X, m["T_curve"], color=COL[p], label=p)
            ax[k].plot(m["x"], m["T_cal"], "o", color=COL[p], ms=5)
        ax[k].set(xlabel="control radius $r_c$ [mm]" if kind == "ECP" else "bisector distance $d$ [mm]",
                  title="ECP (sphere / disc average)" if kind == "ECP" else "point method (TCD PM)",
                  ylim=(1, 3), xlim=(0, 0.5))
    ax[0].set_ylabel("calibration scatter $T$ (V-notch + hole 8)")
    ax[1].legend(fontsize=7)
    _save(fig, OUT, "ecp_compare_calibration")

    tags = [f"{p}, {k}" for p in COL for k in ("ECP", "PM")]
    groups = list(GROUPS)
    fig, axs = plt.subplots(1, 2, figsize=(10.0, 3.8), sharey=True)
    for ax, mode in zip(axs, ("const", "var")):
        for j, tag in enumerate(tags):
            y = np.arange(len(groups)) + (j - (len(tags) - 1) / 2) * 0.12
            g = o["methods"][tag][mode]["groups"]
            p = tag.split(", ")[0]
            ax.errorbar([g[k]["mean"] for k in groups], y, xerr=[g[k]["std"] for k in groups], fmt="o" if "ECP" in tag else "s",
                        mfc=COL[p] if "ECP" in tag else "white", color=COL[p], ms=4, lw=1, label=tag)
        ax.axvline(0, color="#0b0b0b", lw=1)
        for sg in (1, -1):
            ax.axvline(sg * np.log10(3), color="#52514e", lw=0.8, ls="--")
        ax.set(xlabel="$E_N$ = log$_{10}$($N_f/N_{exp}$), mean ± std", xlim=(-1.6, 1.6),
               title="constant length" if mode == "const" else "length = A (N / 1e5)$^B$")
        ax.grid(axis="y", visible=False)
    axs[0].set_yticks(range(len(groups)), groups)
    axs[0].invert_yaxis()
    axs[1].legend(fontsize=7, loc="upper left", bbox_to_anchor=(1.01, 1.0))
    _save(fig, OUT, "ecp_compare_families")


def report(o):
    M = o["methods"]
    f = lambda d: f"{d['mean']:+.2f} / {d['std']:.2f} / {100 * d['in3']:.0f} %"
    rank = sorted(((tag, mode, M[tag][mode]["groups"]["all 19"]["rms"]) for tag in M for mode in ("const", "var")),
                  key=lambda r: r[2])
    best = rank[0]
    pm, ecp = M["FS (length/4), PM"]["const"], M["FS (length/4), ECP"]["const"]
    r = o["swt_repro"]
    L = ["# ECP vs point method, FS vs SWT: our pipeline on the Chiocca & Frendo / Susmel & Taylor En3B data", "",
         "**Summary.**",
         f"* **Best overall: {best[0]}, {'constant' if best[1] == 'const' else 'life-dependent'} length** "
         f"(rms E_N {best[2]:.2f} decades on the 19 prevision series).",
         f"* **Averaging vs point (FS, path-length amplitude, constant length).**",
         f"  * ECP: r_c = {ecp['x']:.3f} mm, E_N {f(ecp['groups']['all 19'])}.",
         f"  * Point method: d = {pm['x']:.3f} mm, E_N {f(pm['groups']['all 19'])}.",
         f"* **SWT is the wrong parameter for this steel, in both forms** (all-19 rms: ECP "
         f"{M['SWT, ECP']['const']['groups']['all 19']['rms']:.2f}, PM {M['SWT, PM']['const']['groups']['all 19']['rms']:.2f}; "
         f"FS ~{M['FS (length/4), ECP']['const']['groups']['all 19']['rms']:.2f}). "
         "It is blind to torsion-dominated cracking, as the paper found.",
         "* **Reproduction of the paper's SWT result with our FE factors (r_c = 0.25 mm):**",
         f"  * calibration T {r['T_cal']:.2f} (paper 1.32);",
         f"  * prevision self-fit T {r['T_prev']:.2f} (paper 4.55);",
         f"  * uniaxial E_N {r['uni']['mean']:+.2f} ± {r['uni']['std']:.2f} (paper +0.30 ± 0.26);",
         f"  * multiaxial E_N {r['multi']['mean']:+.2f} ± {r['multi']['std']:.2f} (paper -0.49 ± 0.65).",
         "",
         "## Calibrated lengths and scatter",
         "| method | x const [mm] | T cal | x(N): A [mm], B | x at 1e4 / 1e6 [mm] | T cal (x(N)) |",
         "|---|---|---|---|---|---|"]
    for tag, d in M.items():
        c, v = d["const"], d["var"]
        L.append(f"| {tag} | {c['x']:.3f} | {c['T_cal']:.3f} | {v['A']:.3f}, {v['B']:+.3f} | "
                 f"{v['x_at']['N1e4']:.3f} / {v['x_at']['N1e6']:.3f} | {v['T_cal']:.3f} |")
    for mode, ttl in (("const", "constant length"), ("var", "life-dependent length x(N)")):
        L += ["", f"## Prediction, {ttl} (mean E_N / std / % within factor 3)",
              "| method | " + " | ".join(GROUPS) + " | rms all |", "|---|" + "---|" * (len(GROUPS) + 1)]
        for tag, d in M.items():
            g = d[mode]["groups"]
            L.append(f"| {tag} | " + " | ".join(f(g[k]) for k in GROUPS) + f" | {g['all 19']['rms']:.2f} |")
    L += ["", "## Paper Table 4 (SWT, r_c = 0.25 mm): per-series self-fit T, ours vs paper",
          "| series | ours | paper |", "|---|---|---|"]
    for n, t in r["per"].items():
        L.append(f"| {n} | {t:.3f} | {PAPER_SWT[n]:.3f} |")
    L += ["", "## Checklist",
          "* **Theory.**",
          "  * FS (extended plane; MCC or path length / 4 shear amplitude) and SWT (Socie: max normal-strain "
          "amplitude plane, σ_n,max on it), both on linear-elastic tensors.",
          "  * ECP averages the tensor over a sphere or disc of radius r_c. PM takes the tensor at distance d "
          "along the notch bisector, at the hot-spot depth for the plates.",
          "* **Physics.**",
          "  * En3B is ductile; FS is the physically appropriate parameter for its shear cracking.",
          "  * Plasticity is not modelled. Calibration uses N_f > 3e4 only, but the prediction set includes "
          "short lives.",
          "* **Numerics.**",
          "  * Fields and meshes as in `ECP_FE.md` (converged at the calibrated radius).",
          "  * Path: 36 samples per cycle.",
          "  * SWT plane: 4000-normal grid plus two cap refinements (to 0.03 deg).",
          "  * x(N): 50 × 37 grid on (A, B). Predicted life is the first crossing of P(x(N)) = c N^m.",
          "  * The V-notch plate is finely meshed only to 0.35 mm.",
          "* **Visualisation.**",
          "  * `ecp_compare_calibration`: calibration scatter vs length for each method.",
          "  * `ecp_compare_families`: prediction error by loading family, constant vs life-dependent length.",
          ""]
    with open(os.path.join(OUT, "ECP_COMPARE.md"), "w") as fh:
        fh.write("\n".join(L) + "\n")


if __name__ == "__main__":
    run()
