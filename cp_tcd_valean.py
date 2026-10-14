"""Blind notch-fatigue prediction and data collapse on Valean et al. (TAFM 131, 2024, 104417): PLA and PLA-CF.

    python cp_tcd_valean.py   -> results/CP_TCD_VALEAN.md, cp_tcd_valean.json, cp_tcd_valean.png

Why this dataset: same printing parameters for smooth and notched (smooth data from the authors' earlier study); individual fatigue results
(15 + 15 per material, R = 0.1) and static notch data for four V opening angles are tabulated.

Model: 3D 1/8 models (thickness 4 mm, mid-plane symmetry), NGSolve P2, unit NET stress (net width 10 mm):
double V-notch (depth 5, opening 2a, root radius rho; not reported -> 0.1 mm, sensitivity 0.05 / 0.2) and the
smooth hourglass (R 45, minimum width 10). E, nu = 0.35 (nu not reported; it barely affects the stress field).
Effective stress: maximum principal stress at the mid-plane on the notch bisector at distance L/2 (TCD point method).
At a single load ratio every normal-stress critical-plane parameter (max normal, SWT, Walker) reduces to this, so the
test here is the notch transfer, not the criterion.

Protocol
  STATIC CALIBRATION   (sigma_0, L_s) from the 4 notch angles x 5 static tests (least squares in log load).
  BLIND PREDICTION     notched fatigue lives from the smooth fatigue curve with L = L_s (no notched fatigue data used).
  COLLAPSE             smooth + notched fatigue in effective stress, one Basquin curve, scatter index
                       T_sigma = 10^(2 * 1.645 * s_log10(sigma)) (5-95 %), for L = L_s (blind) and for the L that
                       minimises T (calibrated, the paper's SED procedure; their T: PLA 1.24, PLA-CF 1.29).
"""
import json
import os
import time

import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

import valean_data as V
from plots import _save

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "results")
NU = 0.35
XR = V.W / 2 - V.DEPTH                                                     # notch root / minimum section at x = 5
LH = V.LEN / 2
ANGLES = (30, 60, 90, 120)
RHOS = (0.05, 0.1, 0.2)
RHO0 = 0.1


def face(kind, alpha=90.0, rho=RHO0):
    """Quarter plate x in [0, W/2], y in [0, LEN/2]; feature at x = W/2 centred on y = 0, root/minimum at x = XR."""
    from netgen.occ import MoveTo, Circle
    f = MoveTo(0, 0).LineTo(V.W / 2, 0).LineTo(V.W / 2, LH).LineTo(0, LH).Close().Face()
    if kind == "smooth":
        return f - Circle((XR + V.R_HOURGLASS, 0), V.R_HOURGLASS).Face()
    a = np.radians(alpha / 2)
    xc = XR + rho
    xt, yt = xc - rho * np.sin(a), rho * np.cos(a)
    xe = V.W / 2 + 1.0
    wedge = MoveTo(xc, 0).LineTo(xe, 0).LineTo(xe, yt + (xe - xt) * np.tan(a)).LineTo(xt, yt).Close().Face()
    return f - wedge - Circle((xc, 0), rho).Face()


def solid_mesh(kind, alpha, rho, scale):
    from netgen.occ import OCCGeometry, Z, X, Y
    s = face(kind, alpha, rho).Extrude(V.T / 2 * Z)
    s.faces.Min(X).name = "xsym"
    s.faces.Min(Y).name = "ysym"
    s.faces.Min(Z).name = "zsym"
    s.faces.Max(Y).name = "load"
    for fc in s.faces:
        if fc.name in ("xsym", "ysym", "zsym", "load"):
            continue
        c = fc.center
        if kind == "notch" and np.hypot(c.x - (XR + rho), c.y) < 3 * rho:
            fc.maxh = scale * rho / 4                                      # root fillet
        elif c.x > XR - 1e-9 and c.y < 6:
            fc.maxh = scale * (0.4 if kind == "notch" else 0.6)
    for e in s.edges:
        c = e.center
        if abs(c.y) < 1e-9 and abs(c.z) < 1e-9:
            e.maxh = scale * (0.04 if kind == "notch" else 0.15)          # bisector (TCD path)
    s.faces.Min(Y).maxh = scale * 0.3
    return OCCGeometry(s).GenerateMesh(maxh=3.0, grading=0.25)


def bisector(kind, alpha, rho, scale, rr, E):
    from ngsolve import Mesh, VectorH1, BilinearForm, LinearForm, GridFunction, Sym, Grad, Id, Trace, InnerProduct, \
        ds, dx, TaskManager, SetNumThreads
    SetNumThreads(4)
    t = time.perf_counter()
    mesh = Mesh(solid_mesh(kind, alpha, rho, scale))
    mu, lam = E / 2 / (1 + NU), E * NU / (1 + NU) / (1 - 2 * NU)
    fes = VectorH1(mesh, order=2, dirichletx="xsym", dirichlety="ysym", dirichletz="zsym")
    u, v = fes.TnT()
    sig = lambda w: 2 * mu * Sym(Grad(w)) + lam * Trace(Sym(Grad(w))) * Id(3)
    with TaskManager():
        a = BilinearForm(InnerProduct(sig(u), Sym(Grad(v))) * dx).Assemble()
        f = LinearForm((V.NET / V.W) * v[1] * ds("load")).Assemble()       # unit NET stress
        g = GridFunction(fes)
        g.vec.data = a.mat.Inverse(fes.FreeDofs(), inverse="sparsecholesky") * f.vec
    S = np.array([sig(g)(mesh(XR - x, 1e-7, 1e-7)) for x in rr]).reshape(len(rr), 3, 3)
    p = np.linalg.eigvalsh(S)[:, 2]
    gx, gw = np.polynomial.legendre.leggauss(32)                          # equilibrium check on y = 0 face
    F = sum(wi * wz * XR / 2 * V.T / 4 * sig(g)[1, 1](mesh(XR * (xi + 1) / 2 - 1e-9, 1e-7, V.T / 4 * (zi + 1)))
            for xi, wi in zip(gx, gw) for zi, wz in zip(*np.polynomial.legendre.leggauss(6)))
    return dict(p=p.tolist(), syy=S[:, 1, 1].tolist(), ndof=int(sum(fes.FreeDofs())), t=time.perf_counter() - t,
                force=float(F / (XR * V.T / 2)))


# ------------------------------------------------------------------------------------------------ TCD analysis
def interp(rr, p, r):
    return np.interp(r, rr, p)


def static_fit(rr, P, mat):
    """sigma_0 and L_s from sigma_net,crit * p_angle(L/2) = sigma_0 over all static tests (log least squares)."""
    best = None
    for L in np.geomspace(0.02, 4.0, 600):
        x = np.array([np.log(Fc / (V.NET * V.T) * interp(rr, P[a], L / 2)) for a in ANGLES for Fc in V.STATIC[mat][a]])
        s0, res = np.exp(x.mean()), np.sqrt(np.mean((x - x.mean()) ** 2))
        if best is None or res < best[2]:
            best = (L, s0, res)
    L, s0, res = best
    pred = {a: float(s0 / interp(rr, P[a], L / 2) * V.NET * V.T) for a in ANGLES}
    return dict(L=float(L), sigma0=float(s0), rms_log=float(res), Fpred=pred,
                Fexp={a: float(np.mean(V.STATIC[mat][a])) for a in ANGLES})


def fatigue_arrays(mat, which):
    d = np.array([(r[1], r[2], r[3]) for r in V.FATIGUE[mat][which]], float)
    return d[:, 0], d[:, 1], d[:, 2].astype(bool)


def basquin(S, N):
    """log10 N = a - k log10 S (least squares on failures); returns a, k and scatter index T_sigma (5-95 %)."""
    k, a = np.polyfit(np.log10(S), np.log10(N), 1)
    k = -k
    s_logS = np.std(np.log10(N) - (a - k * np.log10(S)), ddof=2) / k
    return float(a), float(k), float(10 ** (2 * 1.645 * s_logS))


def collapse(rr, p_s, p_n, mat, L):
    Ss, Ns, ros = fatigue_arrays(mat, "smooth")
    Sn, Nn, ron = fatigue_arrays(mat, "notched")
    es, en = Ss * interp(rr, p_s, L / 2), Sn * interp(rr, p_n, L / 2)
    S = np.concatenate([es[~ros], en[~ron]])
    N = np.concatenate([Ns[~ros], Nn[~ron]])
    a, k, T = basquin(S, N)
    return dict(L=float(L), a=a, k=k, T=T, eff_smooth=es.tolist(), eff_notched=en.tolist())


def blind(rr, p_s, p_n, mat, L):
    """Smooth curve in effective stress (L), notched lives predicted from it."""
    Ss, Ns, ros = fatigue_arrays(mat, "smooth")
    Sn, Nn, ron = fatigue_arrays(mat, "notched")
    a, k, Ts = basquin(Ss[~ros] * interp(rr, p_s, L / 2), Ns[~ros])
    Np = 10 ** (a - k * np.log10(Sn * interp(rr, p_n, L / 2)))
    lf = np.log10(Np[~ron] / Nn[~ron])
    return dict(a=a, k=k, T_smooth=Ts, N_pred=Np.tolist(), rms_log=float(np.sqrt(np.mean(lf ** 2))),
                mean_log=float(lf.mean()), in3=float(np.mean(np.abs(lf) < np.log10(3))))


def split_validation(rr, p_s, p_n, mat):
    """L calibrated on smooth data + notched failures with N below the notched median; the notched lives above the
    median are then predicted from the smooth curve (they are never used)."""
    Ss, Ns, ros = fatigue_arrays(mat, "smooth")
    Sn, Nn, ron = fatigue_arrays(mat, "notched")
    med = np.median(Nn[~ron])
    cal = (~ron) & (Nn < med)
    test = (~ron) & (Nn >= med)

    def T_of(L):
        S = np.concatenate([Ss[~ros] * interp(rr, p_s, L / 2), Sn[cal] * interp(rr, p_n, L / 2)])
        return basquin(S, np.concatenate([Ns[~ros], Nn[cal]]))[2]
    L = float(min(np.geomspace(0.05, 4.0, 400), key=T_of))
    a, k, _ = basquin(Ss[~ros] * interp(rr, p_s, L / 2), Ns[~ros])
    lf = np.log10(10 ** (a - k * np.log10(Sn[test] * interp(rr, p_n, L / 2))) / Nn[test])
    return dict(L=L, n_cal=int(cal.sum()), n_test=int(test.sum()), rms_log=float(np.sqrt(np.mean(lf ** 2))),
                mean_log=float(lf.mean()), in3=float(np.mean(np.abs(lf) < np.log10(3))), N_split=float(med))


def main():
    rr = np.concatenate([[0.0], np.geomspace(1e-3, 0.05, 12), np.linspace(0.06, 4.9, 160)])
    out = dict(rr=rr.tolist(), fe={})
    for mat in V.TENSILE:
        E = V.TENSILE[mat]["E"]
        fe = {}
        for rho in RHOS:
            for a in ANGLES:
                d = bisector("notch", a, rho, 1.0, rr, E)
                fe[f"notch|{a}|{rho}"] = d
                print(mat, "notch", a, rho, "ndof", d["ndof"], "t %.1f" % d["t"], "Kt(sigma1) %.2f" % d["p"][0],
                      "force/net %.4f" % d["force"], flush=True)
        d = bisector("notch", 90, RHO0, 0.5, rr, E)
        fe[f"notch|90|{RHO0}|fine"] = d
        print(mat, "notch 90 fine ndof", d["ndof"], "Kt %.3f" % d["p"][0], flush=True)
        fe["smooth"] = bisector("smooth", 0, 0, 1.0, rr, E)
        print(mat, "smooth Kt %.3f force %.4f" % (fe["smooth"]["p"][0], fe["smooth"]["force"]), flush=True)
        out["fe"][mat] = fe
    res = {}
    for mat in V.TENSILE:
        fe = out["fe"][mat]
        p_s = np.array(fe["smooth"]["p"])
        r = {}
        for rho in RHOS:
            P = {a: np.array(fe[f"notch|{a}|{rho}"]["p"]) for a in ANGLES}
            st = static_fit(rr, P, mat)
            bl = blind(rr, p_s, P[90], mat, st["L"])
            col_blind = collapse(rr, p_s, P[90], mat, st["L"])
            Ls = np.geomspace(0.05, 4.0, 400)
            Topt = min(Ls, key=lambda L: collapse(rr, p_s, P[90], mat, L)["T"])
            col_opt = collapse(rr, p_s, P[90], mat, Topt)
            sv = split_validation(rr, p_s, P[90], mat)
            r[str(rho)] = dict(static=st, blind=bl, collapse_blind=col_blind, collapse_opt=col_opt, split=sv)
            print(mat, rho, "split: L %.2f, predicted long-life half rms log %.2f mean %+.2f in3 %.0f%%" % (
                sv["L"], sv["rms_log"], sv["mean_log"], 100 * sv["in3"]), flush=True)
            print(mat, rho, "L_s %.3f sigma0 %.1f (UTS %.1f) static rms %.3f | blind rms log %.2f mean %+.2f | "
                  "T blind %.2f  T opt %.2f (L %.2f)" % (st["L"], st["sigma0"], V.TENSILE[mat]["UTS"], st["rms_log"],
                                                        bl["rms_log"], bl["mean_log"], col_blind["T"], col_opt["T"],
                                                        col_opt["L"]), flush=True)
        Sn, Nn, ron = fatigue_arrays(mat, "notched")
        Ss, Ns, ros = fatigue_arrays(mat, "smooth")
        r["nominal"] = dict(smooth=basquin(Ss[~ros], Ns[~ros]), notched=basquin(Sn[~ron], Nn[~ron]),
                            combined=basquin(np.concatenate([Ss[~ros], Sn[~ron]]), np.concatenate([Ns[~ros], Nn[~ron]])))
        pf, pc = np.array(fe[f"notch|90|{RHO0}|fine"]["p"]), np.array(fe[f"notch|90|{RHO0}"]["p"])
        L = r[str(RHO0)]["collapse_opt"]["L"]
        r["mesh_sens"] = float(abs(interp(rr, pc, L / 2) / interp(rr, pf, L / 2) - 1))
        res[mat] = r
    out["res"] = res
    json.dump(out, open(os.path.join(OUT, "cp_tcd_valean.json"), "w"), indent=1)
    return out


def t_convention(res):
    """Our T = 10^(2*1.645*s); the paper's nominal T values are systematically larger. Exponent f with
    log T_paper = f log T_ours, fitted on the paper's 4 nominal curves (smooth / notched x 2 materials)."""
    ours, theirs = [], []
    for m in res:
        for g in ("smooth", "notched"):
            ours.append(np.log(res[m]["nominal"][g][2]))
            theirs.append(np.log(V.PAPER[m]["T_nom"][g]))
    ours, theirs = np.array(ours), np.array(theirs)
    return float(np.sum(ours * theirs) / np.sum(ours ** 2)), (theirs / ours).tolist()


def report(out=None):
    out = out or json.load(open(os.path.join(OUT, "cp_tcd_valean.json")))
    res, fe = out["res"], out["fe"]
    L = ["# Blind notch-fatigue prediction and data collapse: Valean et al. (TAFM 131, 2024) PLA and PLA-CF", "",
         "Data: Valean, Foti, Berto, Marsavina, Theor. Appl. Fract. Mech. 131 (2024) 104417 (CC BY), Tables 2-4: static "
         "V-notch critical loads (2a = 30/60/90/120 deg, 5 each) and individual fatigue results (15 smooth + 15 notched "
         "2a = 90 deg, R = 0.1; same printing parameters, smooth fatigue data from the authors' earlier study [39]). Model: 3D 1/8, NGSolve P2, max principal stress at the mid-plane on "
         "the bisector; root radius not reported -> rho = 0.1 mm (sensitivity 0.05 / 0.2 mm). Static calibration "
         "(sigma_0, L_s) on the 4 angles; notched fatigue then predicted blind from the smooth curve with L_s.", "",
         "## 1. Checks", "", "| material | geometry | DOF | Kt (sigma1, net) | ligament force / net load |", "|---|---|---|---|---|"]
    for mat in fe:
        for a in ANGLES:
            d = fe[mat][f"notch|{a}|{RHO0}"]
            L.append("| %s | V %d deg, rho %.2f | %d | %.2f | %.4f |" % (mat, a, RHO0, d["ndof"], d["p"][0], d["force"]))
        d = fe[mat]["smooth"]
        L.append("| %s | smooth hourglass R 45 | %d | %.3f | %.4f |" % (mat, d["ndof"], d["p"][0], d["force"]))
    L += ["", "Mesh sensitivity of sigma1 at L/2 (fatigue L, V 90, rho 0.1, coarse vs fine): " + ", ".join(
        "%s %.1e" % (m, res[m]["mesh_sens"]) for m in res) + ".", "",
          "## 2. Static TCD calibration (4 angles x 5 tests)", "",
          "| material | rho [mm] | L_s [mm] | sigma_0 [MPa] (UTS) | rms log load | F_pred / F_exp: 30 / 60 / 90 / 120 deg |",
          "|---|---|---|---|---|---|"]
    for mat in res:
        for rho in RHOS:
            st = res[mat][str(rho)]["static"]
            L.append("| %s | %.2f | %.2f | %.1f (%.1f) | %.3f | %s |" % (
                mat, rho, st["L"], st["sigma0"], V.TENSILE[mat]["UTS"], st["rms_log"],
                " / ".join("%.2f" % (st["Fpred"][str(a)] / st["Fexp"][str(a)]) if isinstance(list(st["Fpred"])[0], str)
                           else "%.2f" % (st["Fpred"][a] / st["Fexp"][a]) for a in ANGLES)))
    L += ["", "Paper's SED control radius: PLA R0 = 0.81 mm (static), 0.91-0.97 mm (fatigue); PLA-CF 0.30 / 0.52-0.69 mm.", "",
          "## 3. Fatigue: collapse, split-sample prediction, static-L blind prediction", "",
          "| material | rho | T_sigma nominal (smooth+notched) | T_sigma TCD, L optimised (L) | **same, paper's convention** | paper SED T_sigma | "
          "T_sigma smooth alone / notched alone | split test: long-life notched predicted, rms log10 / mean / within x3 (L) | "
          "static-L blind: rms log10 / mean / within x3 | T_sigma with static L |", "|---|---|---|---|---|---|---|---|---|---|"]
    conv = t_convention(res)
    for mat in res:
        nom = res[mat]["nominal"]
        for rho in RHOS:
            q = res[mat][str(rho)]
            sv = q["split"]
            L.append("| %s | %.2f | %.2f | %.2f (%.2f mm) | **%.2f** | %.2f | %.2f / %.2f | %.2f / %+.2f / %.0f %% (%.2f mm) | "
                     "%.2f / %+.2f / %.0f %% | %.2f |" % (
                         mat, rho, nom["combined"][2], q["collapse_opt"]["T"], q["collapse_opt"]["L"],
                         q["collapse_opt"]["T"] ** conv[0], V.PAPER[mat]["T_SED"],
                         nom["smooth"][2], nom["notched"][2], sv["rms_log"], sv["mean_log"], 100 * sv["in3"], sv["L"],
                         q["blind"]["rms_log"], q["blind"]["mean_log"], 100 * q["blind"]["in3"], q["collapse_blind"]["T"]))
    L += ["", "**Scatter-index convention.** Our T = 10^(2 x 1.645 s); on the paper's own nominal curves it gives "
          "values systematically below the paper's (ratio of logs per curve: " + ", ".join("%.2f" % x for x in conv[1]) +
          "; least-squares exponent %.2f), i.e. the paper uses a wider quantile or a different estimator. Converted to "
          "the paper's convention, one fatigue L collapses smooth + notched about **as well as, or slightly worse than, "
          "the paper's SED** (PLA ~1.29 vs 1.24, PLA-CF ~1.36 vs 1.29) - not better." % conv[0], "",
          "**What the collapse is.** With a single fatigue notch geometry, the optimised L acts as a fitted fatigue notch "
          "factor (T_opt is identical for every assumed rho): it shows that one constant reconciles smooth and notched "
          "data, not that TCD transfers between notch geometries. The split test splits by stress level within the same "
          "geometry: a same-geometry extrapolation check (mostly that smooth and notched slopes are similar, 5.43 vs 5.10 "
          "for PLA), not a prediction for a new notch.", "",
          "**Static calibration works; static -> fatigue transfer of L does not, as expected.** The linear-elastic point "
          "method with (sigma_0, L_s) fitted on the 4 angles reproduces the static loads to ~7 % rms, including the minimum "
          "at 60 deg (comparable to the paper's SED, +-7 %); sigma_0 > UTS is normal for linear-elastic TCD on polymers. "
          "But L_s is not identifiable with the root radius unknown (0.10-0.30 mm for rho 0.2-0.05) and L is life-dependent "
          "(the paper's own R0 grows from static to fatigue: 0.81 -> 0.97 mm PLA, 0.30 -> 0.69 mm PLA-CF), so using L_s "
          "for high-cycle fatigue is not a valid TCD prediction - TCD does not claim it."]
    L += ["", "Nominal Basquin fits (ours vs paper, log10 N = a - k log10 dS): " + "; ".join(
        "%s smooth %.2f / %.2f (paper %.2f / %.2f), notched %.2f / %.2f (paper %.2f / %.2f)" % (
            m, res[m]["nominal"]["smooth"][0], res[m]["nominal"]["smooth"][1], *V.PAPER[m]["smooth"],
            res[m]["nominal"]["notched"][0], res[m]["nominal"]["notched"][1], *V.PAPER[m]["notched"]) for m in res) + ".",
          "", "![results](cp_tcd_valean.png)"]
    open(os.path.join(OUT, "CP_TCD_VALEAN.md"), "w").write("\n".join(L) + "\n")
    fig(out)


def fig(out):
    rr = np.array(out["rr"])
    res, fe = out["res"], out["fe"]
    f, ax = plt.subplots(2, 4, figsize=(19, 8.8), layout="constrained")
    cmap = plt.get_cmap("viridis")
    for i, mat in enumerate(res):
        q = res[mat][str(RHO0)]
        Ls = q["static"]["L"]
        a = ax[i, 0]
        for j, ang in enumerate(ANGLES):
            a.plot(rr, fe[mat][f"notch|{ang}|{RHO0}"]["p"], color=cmap(j / 3.5), lw=1.5, label=f"V {ang}°")
        a.plot(rr, fe[mat]["smooth"]["p"], "k--", lw=1.2, label="smooth hourglass")
        Lo = q["collapse_opt"]["L"]
        a.axvline(Lo / 2, color="0.4", ls=":", lw=1)
        a.text(Lo / 2, 4.6, " L/2 (fatigue)", fontsize=8, color="0.3")
        a.set(xlim=(0, 3), ylim=(0, 5), xlabel="distance from root along bisector r [mm]",
              ylabel=r"$\sigma_1/\sigma_{net}$ (mid-plane)", title=f"({'ae'[i]}) {mat}: bisector stress, ρ = {RHO0} mm")
        a.legend(fontsize=7)
        a.grid(alpha=0.3)
        a = ax[i, 1]
        for j, ang in enumerate(ANGLES):
            a.plot([ang] * 5, V.STATIC[mat][ang], "o", color="0.6", ms=4, mfc="none", label="tests" if j == 0 else None)
        st = q["static"]
        a.plot(ANGLES, [st["Fpred"][str(x)] if str(x) in st["Fpred"] else st["Fpred"][x] for x in ANGLES], "rs-",
               label=f"TCD-PM fit: L_s = {Ls:.2f} mm, σ₀ = {st['sigma0']:.0f} MPa")
        a.set(xlabel="notch opening angle 2α [deg]", ylabel="critical load [N]", xticks=ANGLES,
              title=f"({'bf'[i]}) static calibration (rms log {st['rms_log']:.3f})")
        a.legend(fontsize=7)
        a.grid(alpha=0.3)
        Ss, Ns, ros = fatigue_arrays(mat, "smooth")
        Sn, Nn, ron = fatigue_arrays(mat, "notched")
        a = ax[i, 2]
        for S, N, ro, mk, lab in ((Ss, Ns, ros, "o", "smooth"), (Sn, Nn, ron, "s", "notched V 90°")):
            a.loglog(N[~ro], S[~ro], mk, color="#2e86c1" if mk == "o" else "#c0392b", ms=5, label=lab)
            a.loglog(N[ro], S[ro], mk, color="#2e86c1" if mk == "o" else "#c0392b", ms=5, mfc="none")
        T = res[mat]["nominal"]["combined"][2]
        a.set(xlabel="cycles to failure N", ylabel="Δσ nominal (net) [MPa]", ylim=(2, 60),
              title=f"({'cg'[i]}) nominal: two curves, combined Tσ = {T:.2f}")
        a.legend(fontsize=7)
        a.grid(alpha=0.3, which="both")
        a = ax[i, 3]
        cb = q["collapse_opt"]
        es, en = np.array(cb["eff_smooth"]), np.array(cb["eff_notched"])
        a.loglog(Ns[~ros], es[~ros], "o", color="#2e86c1", ms=5, label="smooth")
        a.loglog(Ns[ros], es[ros], "o", color="#2e86c1", ms=5, mfc="none")
        a.loglog(Nn[~ron], en[~ron], "s", color="#c0392b", ms=5, label="notched V 90\u00b0")
        a.loglog(Nn[ron], en[ron], "s", color="#c0392b", ms=5, mfc="none")
        Ng = np.geomspace(3e2, 3e6, 50)
        Sg = 10 ** ((cb["a"] - np.log10(Ng)) / cb["k"])
        Tb = cb["T"]
        a.loglog(Ng, Sg, "k-", lw=1)
        a.fill_between(Ng, Sg / np.sqrt(Tb), Sg * np.sqrt(Tb), color="0.9", label=f"5-95 %: Tσ = {Tb:.2f}")
        a.set(xlabel="cycles to failure N", ylabel=r"$\Delta\sigma_{eff}=\Delta\sigma_1(L/2)$ [MPa]", ylim=(2, 60),
              title=f"({'dh'[i]}) TCD collapse, one L = {cb['L']:.2f} mm (paper SED: T\u03c3 = {V.PAPER[mat]['T_SED']:.2f})")
        a.legend(fontsize=7)
        a.grid(alpha=0.3, which="both")
        for a in ax[i]:
            a.title.set_fontsize(9)
    _save(f, OUT, "cp_tcd_valean")


if __name__ == "__main__":
    report(main())
