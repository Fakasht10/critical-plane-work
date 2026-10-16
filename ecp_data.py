"""Loader for the ECP dataset of Chiocca & Frendo, Int. J. Fatigue 189 (2024) 108565 (github.com/achiocca1/ECP).

Material: En3B low-carbon steel (paper Table 1): E = 200 GPa, nu = 0.3, sigma_y = 606.2 MPa, Fatemi-Socie k = 0.4.
Fatigue data: Susmel & Taylor, FFEMS 30 (2007) 567 (flat notched plates, axial, R = -1 / 0.1) and FFEMS 31 (2008)
1047 (cylindrical V-notched bars, tension / torsion / in-phase and 90 deg out-of-phase, R = -1 / 0).

The ECP repository is GPL-3; its files are not copied here. `fetch()` clones it at a pinned commit into data/ecp
(git-ignored). This module only parses those files:

  1_<name>.txt     uniaxial / torsion: 2 rows = degree-6 polynomials in r_c [mm] (MATLAB order, highest power first)
                   of the sphere-averaged (3D) or disc-averaged (2D) Delta gamma_max and sigma_n,max per unit load;
                   multiaxial: one row per specimen = polynomial of 2 FS(r_c) at that specimen's loads.
  2_Exp_<name>.txt N_f, nominal stress amplitude [MPa], load [N or N mm] (multiaxial: N_f, sigma_a, tau_a).
  3_Corr_<name>.txt one factor per specimen: MPHP / MCC ratio of Delta gamma for the 90 deg out-of-phase series.

FS(r_c) = 0.5 (F Dg + F^2 (k / sigma_y) Dg sn)  (paper Eq. 7, F = load), or 0.5 * polynomial for multiaxial rows.
"""
import os
import re
import subprocess

import numpy as np

URL, COMMIT = "https://github.com/achiocca1/ECP", "9d5115258d5a13b926788a0c636a417d94e85efd"
HERE = os.path.dirname(os.path.abspath(__file__))
DIR = os.environ.get("ECP_DIR", os.path.join(HERE, "data", "ecp"))
K_FS, SY, E, NU = 0.4, 606.2, 200e3, 0.3


def fetch(path=DIR):
    if not os.path.isdir(os.path.join(path, "PREVISION")):
        os.makedirs(os.path.dirname(path), exist_ok=True)
        subprocess.run(["git", "clone", "-q", URL, path], check=True)
    subprocess.run(["git", "-C", path, "checkout", "-q", COMMIT], check=True)
    return path


def _parse(name):
    """'Axial_V_notch_R_0.1' -> dict(kind, geom, phase, R)."""
    m = re.match(r"(Axial|Torsion|Multi)_(.+?)(?:_Ph_(\d+))?_R_(-?[\d.]+)$", name)
    kind, geom, ph, R = m.groups()
    return dict(kind=kind, geom=geom, phase=None if ph is None else int(ph), R=float(R))


LABEL = {"V_notch": "V-notch (flat)", "U_notch": "U-notch (flat)", "Hole_8": "Hole 8 mm", "Hole_3_5": "Hole 3.5 mm",
         "V_notch_02": "V-notch r=0.2", "V_notch_12": "V-notch r=1.25", "V_notch_4": "V-notch r=4"}


def load(path=None, k=K_FS, sy=SY, corrector=True):
    """List of series dicts: name, role ('calibration'|'prevision'), kind, geom, phase, R, Nf, S, load, label,
    fs(rc) -> (n_specimens, len(rc)) ECP Fatemi-Socie factor, and for uniaxial/torsion dg(rc), sn(rc) per unit load."""
    path = path or fetch()
    out = []
    for role, sub in (("calibration", "CALIBRATION"), ("prevision", "PREVISION")):
        d = os.path.join(path, sub)
        for fn in sorted(os.listdir(d)):
            if not fn.startswith("2_Exp_"):
                continue
            name = fn[6:-4]
            exp = np.loadtxt(os.path.join(d, fn), ndmin=2)
            P = np.loadtxt(os.path.join(d, "1_" + name + ".txt"), ndmin=2)
            cf = os.path.join(path, "CORRECTOR", "3_Corr_" + name + ".txt")
            corr = np.loadtxt(cf, ndmin=1) if (corrector and os.path.isfile(cf)) else np.ones(len(exp))
            s = dict(name=name, role=role, **_parse(name), Nf=exp[:, 0], S=exp[:, 1], load=exp[:, 2], P=P,
                     corr=corr, multi=P.shape[0] != 2)
            if s["multi"]:
                assert P.shape[0] == len(exp), name
                s["tau"] = exp[:, 2]
            s["label"] = LABEL[s["geom"]] + (f", {s['phase']} deg" if s["phase"] is not None else "") + \
                (", torsion" if s["kind"] == "Torsion" else "") + f", R={s['R']:g}"
            s["fs"] = _fs_fun(s, k, sy)
            out.append(s)
    return out


def _fs_fun(s, k, sy):
    P, F, corr = s["P"], s["load"], s["corr"]
    if s["multi"]:
        return lambda rc: corr[:, None] * 0.5 * np.array([np.polyval(p, rc) for p in P])
    dg = lambda rc: np.polyval(P[0], rc)
    sn = lambda rc: np.polyval(P[1], rc)
    s["dg"], s["sn"] = dg, sn
    return lambda rc: corr[:, None] * 0.5 * (F[:, None] * dg(rc)[None] + F[:, None] ** 2 * (k / sy) * (dg(rc) * sn(rc))[None])


def reference_results(path=None):
    """The repository's own outputs (MATLAB run, r_c = 0.2): T, c10, c50, c90, m."""
    path = path or fetch()
    r = os.path.join(path, "RESULTS")
    rd = lambda f: np.loadtxt(os.path.join(r, f))
    per = {fn[len("FatigueScatter_FS_1_"):-4]: rd(fn) for fn in os.listdir(r) if fn.startswith("FatigueScatter_FS_1_")}
    return dict(calib=rd("CalibScatter_FS.txt"), total=rd("TotalFatigueScatter_FS.txt"), per=per,
                residuals=rd("FS_Residuals.txt"), radius_vs_nf=rd("radiusvsNf_FS.txt"),
                fs={fn[5:-4]: rd(fn) for fn in os.listdir(r) if fn.startswith("FS_1_")})
