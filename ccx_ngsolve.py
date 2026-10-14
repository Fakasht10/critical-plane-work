"""Toolchain check: NGSolve P2 vs CalculiX C3D10 on the SAME mesh -> nodal stress -> critical-plane search (cplane).

    python ccx_ngsolve.py   -> results/CCX_NGSOLVE.md, ccx_ngsolve.json, ccx_ngsolve.png

Geometry: Ezeh & Susmel (2020) intermediate notch (U, r = 1 mm; W 25, net 6, t 3 mm), 1/8 model (symmetry x = 0,
y = 0 ligament, z = 0 mid-plane), unit gross traction on y = L/2. PLA E = 3476 MPa, nu = 0.33 (ezeh_geometry.py).
Same straight-sided netgen tetrahedra for both solvers; mid-edge nodes added on straight edges, so NGSolve VectorH1
order 2 and ccx C3D10 span the same space and integrate the same stiffness exactly.  Stress recovery is the same too:
the P2 stress is P1 per element, ccx extrapolates its 4 Gauss points (exact for P1) and averages over the elements at
a node; we average NGSolve's element-vertex stresses the same way (unweighted).
Checks: displacement and stress agreement (frd: 5 significant digits), Kt_net mesh convergence (3 levels), critical-
plane damage on the free-surface nodes (Findley k = 0.3, MWCM, SWT-stress; R = -1 and R = 0) from both stress fields,
and wall-clock per stage.
"""
import json
import os
import tempfile
import time

import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

import ezeh_geometry as M
import ccx_io
import cplane as cp
from plots import _save

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "results")
E, NU = M.E, M.NU
R_NOTCH = 1.0
LEVELS = (0.4, 0.2, 0.1)                                                # notch-surface mesh size [mm]
SIG_A = 10.0 * (M.W - 2 * M.C) / M.W                                    # gross amplitude for 10 MPa net amplitude


def solid_mesh(h_notch, maxh=3.0):
    from netgen.occ import OCCGeometry, Z, X, Y
    f, xr = M.plate_face("U", R_NOTCH)
    s = f.Extrude(M.T / 2 * Z)
    s.faces.Min(X).name = "xsym"
    s.faces.Min(Y).name = "ysym"
    s.faces.Min(Z).name = "zsym"
    s.faces.Max(Y).name = "load"
    for fc in s.faces:
        c = fc.center
        if c.x > xr - 1e-9 and c.y < 3 * R_NOTCH + 1e-9 and fc.name not in ("ysym", "zsym"):
            fc.maxh = h_notch
    s.faces.Min(Y).maxh = 2 * h_notch
    ng = OCCGeometry(s).GenerateMesh(maxh=maxh, grading=0.3)
    return ng, xr


def mesh_arrays(ng):
    X = np.array([p.p for p in ng.Points()])
    T = np.array([[v.nr - 1 for v in el.vertices] for el in ng.Elements3D()])
    bnd = {}
    names = ng.GetRegionNames(dim=2)
    for el in ng.Elements2D():
        bnd.setdefault(names[el.index - 1], []).append([v.nr - 1 for v in el.vertices])
    return X, T, {k: np.array(v) for k, v in bnd.items()}


def solve_ngsolve(ng):
    from ngsolve import Mesh, VectorH1, BilinearForm, LinearForm, GridFunction, Sym, Grad, Id, Trace, InnerProduct, \
        ds, dx, IntegrationRule, VOL, TaskManager, SetNumThreads, CF, x, y, z
    SetNumThreads(4)
    mesh = Mesh(ng)
    mu, lam = E / 2 / (1 + NU), E * NU / (1 + NU) / (1 - 2 * NU)
    fes = VectorH1(mesh, order=2, dirichletx="xsym", dirichlety="ysym", dirichletz="zsym")
    u, v = fes.TnT()
    sig = lambda w: 2 * mu * Sym(Grad(w)) + lam * Trace(Sym(Grad(w))) * Id(3)
    t = time.perf_counter()
    with TaskManager():
        a = BilinearForm(InnerProduct(sig(u), Sym(Grad(v))) * dx).Assemble()
        f = LinearForm(v[1] * ds("load")).Assemble()
        gfu = GridFunction(fes)
        gfu.vec.data = a.mat.Inverse(fes.FreeDofs(), inverse="sparsecholesky") * f.vec
    t_solve = time.perf_counter() - t
    t = time.perf_counter()
    ir = IntegrationRule([(1, 0, 0), (0, 1, 0), (0, 0, 1), (0, 0, 0)], [0.25] * 4)  # NGSolve tet vertex order
    pts = mesh.MapToAllElements(ir, VOL)
    S = np.array(sig(gfu)(pts)).reshape(mesh.ne, 4, 9)
    Xe = np.array(CF((x, y, z))(pts)).reshape(mesh.ne, 4, 3)
    Ue = np.array(gfu(pts)).reshape(mesh.ne, 4, 3)
    T = np.array([[v.nr for v in el.vertices] for el in mesh.Elements(VOL)])
    nv = mesh.nv
    Xv = np.array([v.point for v in mesh.vertices])
    assert np.abs(Xe - Xv[T]).max() < 1e-12, "vertex order of the integration rule"
    acc, cnt, U = np.zeros((nv, 9)), np.zeros(nv), np.zeros((nv, 3))
    np.add.at(acc, T.ravel(), S.reshape(-1, 9))
    np.add.at(cnt, T.ravel(), 1)
    U[T.ravel()] = Ue.reshape(-1, 3)
    Sv = acc / cnt[:, None]
    Sv = Sv[:, [0, 4, 8, 1, 5, 6]]                                      # -> [xx, yy, zz, xy, yz, zx]
    return dict(S=Sv, U=U, X=Xv, ndof=sum(fes.FreeDofs()), t_solve=t_solve, t_post=time.perf_counter() - t)


def solve_ccx(X, T, bnd, workdir):
    T4 = ccx_io.orient(X, T)
    X10, T10 = ccx_io.quadratic(X, T4)
    tol = 1e-9
    fix = {1: np.where(np.abs(X10[:, 0]) < tol)[0], 2: np.where(np.abs(X10[:, 1]) < tol)[0],
           3: np.where(np.abs(X10[:, 2]) < tol)[0]}
    load = set(map(tuple, np.sort(bnd["load"], axis=1)))
    faces = [(e, k) for e in range(len(T4)) for k, fc in enumerate(ccx_io.C3D10_FACES)
             if tuple(sorted(T4[e, list(fc)])) in load]
    assert len(faces) == len(load)
    t = time.perf_counter()
    ccx_io.write_c3d10_inp(os.path.join(workdir, "job.inp"), X10, T10, E, NU, fix, faces, pressure=-1.0)
    t_write = time.perf_counter() - t
    t_solve, log = ccx_io.run_ccx(workdir, "job")
    t = time.perf_counter()
    nodes, U, S = ccx_io.read_frd(os.path.join(workdir, "job.frd"))
    t_read = time.perf_counter() - t
    keep = nodes < len(X)                                                # corner nodes = mesh vertices
    Sv, Uv = np.zeros((len(X), 6)), np.zeros((len(X), 3))
    Sv[nodes[keep]], Uv[nodes[keep]] = S[keep], U[keep]
    hdr = [l.split()[1] for l in open(os.path.join(workdir, "job.frd")) if l.startswith(" -5")]
    return dict(S=Sv, U=Uv, ndof=3 * len(X10) - sum(len(v) for v in fix.values()), t_solve=t_solve,
                t_write=t_write, t_read=t_read, frd_components=hdr)


def free_surface_nodes(bnd, nv):
    free = np.zeros(nv, bool)
    for name, F in bnd.items():
        if name not in ("xsym", "ysym", "zsym", "load"):
            free[F.ravel()] = True
    return np.where(free)[0]


def criteria():
    return {"Findley k=0.3": cp.Findley(0.3), "MWCM": cp.MWCM(), "SWT-stress": cp.SWTStress()}


def run_level(h):
    t = time.perf_counter()
    ng, xr = solid_mesh(h)
    t_mesh = time.perf_counter() - t
    ngs = solve_ngsolve(ng)
    X, T, bnd = mesh_arrays(ng)                                          # after NGSolve took the mesh
    assert np.abs(ngs["X"] - X).max() < 1e-12, "vertex numbering differs between NGSolve and the exported mesh"
    with tempfile.TemporaryDirectory() as wd:
        ccx = solve_ccx(X, T, bnd, wd)
    surf = free_surface_nodes(bnd, len(X))
    root = np.argmin(np.where((np.abs(X[:, 1]) < 1e-9), np.abs(X[:, 0] - xr) + np.abs(X[:, 2] - M.T / 2), np.inf))
    rootm = np.argmin(np.where((np.abs(X[:, 1]) < 1e-9), np.abs(X[:, 0] - xr) + np.abs(X[:, 2]), np.inf))
    net = M.W / (M.W - 2 * M.C)
    row = dict(h=h, nv=len(X), ne=len(T), t_mesh=t_mesh, ndof_ngs=int(ngs["ndof"]), ndof_ccx=int(ccx["ndof"]),
               t_ngs=ngs["t_solve"], t_ngs_post=ngs["t_post"], t_ccx=ccx["t_solve"], t_ccx_io=ccx["t_write"] + ccx["t_read"],
               frd_components=ccx["frd_components"],
               Kt_surf_ngs=ngs["S"][root, 1] / net, Kt_surf_ccx=ccx["S"][root, 1] / net,
               Kt_mid_ngs=ngs["S"][rootm, 1] / net, Kt_mid_ccx=ccx["S"][rootm, 1] / net,
               root_xyz=X[root].tolist(), n_surf=len(surf))
    scale = np.abs(ngs["S"]).max()
    row["dS_max_rel"] = float(np.abs(ngs["S"] - ccx["S"]).max() / scale)
    row["dU_max_rel"] = float(np.abs(ngs["U"] - ccx["U"]).max() / np.abs(ngs["U"]).max())
    cpl = {}
    for R in (-1, 0):
        for cname, crit in criteria().items():
            res = {}
            for sname, s in (("ngs", ngs), ("ccx", ccx)):
                Sa = s["S"][surf] * SIG_A
                Sm = Sa * (1.0 if R == 0 else 0.0)
                t = time.perf_counter()
                r = cp.search_prop(Sa, Sm, crit)
                dt = time.perf_counter() - t
                t = time.perf_counter()
                rs = cp.search_prop(Sa, Sm, crit, screen=0.5, chunk=256)
                dts = time.perf_counter() - t
                i = int(np.argmax(r["value"]))
                res[sname] = dict(t=dt, t_screen=dts, n_screen=int(rs["n_searched"]), screen_ok=bool(
                    np.isclose(np.nanmax(rs["value"]), r["value"][i])), vmax=float(r["value"][i]), node=int(surf[i]),
                    normal=r["normal"][i].tolist(), tau_a=float(r["tau_a"][i]), sn_a=float(r["sn_a"][i]),
                    sn_m=float(r["sn_m"][i]), field=r["value"])
            a, b = res["ngs"], res["ccx"]
            cpl[f"R={R}|{cname}"] = dict(
                vmax_ngs=a["vmax"], vmax_ccx=b["vmax"], dv=abs(a["vmax"] - b["vmax"]) / a["vmax"],
                field_drel=float(np.max(np.abs(a["field"] - b["field"])) / a["vmax"]),
                same_node=a["node"] == b["node"],
                dnormal_deg=float(np.degrees(np.arccos(min(1.0, abs(np.dot(a["normal"], b["normal"])))))),
                normal=a["normal"], xyz=X[a["node"]].tolist(), tau_a=a["tau_a"], sn_a=a["sn_a"], sn_m=a["sn_m"],
                t=a["t"], t_screen=a["t_screen"], n_screen=a["n_screen"], screen_ok=a["screen_ok"] and b["screen_ok"])
            print(h, R, cname, {k: v for k, v in cpl[f"R={R}|{cname}"].items() if k not in ("normal", "xyz")}, flush=True)
    row["cplane"] = cpl
    if h == LEVELS[-1]:
        full = np.full(len(X), np.nan)
        full[surf] = res["ngs"]["field"]                                 # last loop item: R = 0, SWT-stress
        tri = np.vstack([F for k, F in bnd.items() if k not in ("xsym", "ysym", "zsym", "load")])
        row["_plot"] = dict(X=X.tolist(), F=full.tolist(), tri=tri.tolist(), xr=xr)
    print({k: v for k, v in row.items() if k not in ("cplane", "_plot")}, flush=True)
    return row


def kt_2d(order):
    """Independent anchor: 2D plane stress on the curved geometry (notch edge h 0.03 mm), order-p elements."""
    from ngsolve import Mesh, VectorH1, BilinearForm, LinearForm, GridFunction, Sym, Grad, Id, Trace, InnerProduct, \
        ds, dx
    from netgen.occ import OCCGeometry
    f, xr = M.plate_face("U", R_NOTCH)
    for e in f.edges:
        c = e.center
        if abs(c.x) < 1e-9:
            e.name = "xsym"
        if abs(c.y) < 1e-9:
            e.name = "ysym"
        if abs(c.y - M.L / 2) < 1e-9:
            e.name = "load"
        if c.x > xr - 1e-9 and c.y < 3:
            e.maxh = 0.03
    mesh = Mesh(OCCGeometry(f, dim=2).GenerateMesh(maxh=2.0, grading=0.3))
    mesh.Curve(order)
    fes = VectorH1(mesh, order=order, dirichletx="xsym", dirichlety="ysym")
    u, v = fes.TnT()
    mu, lam = E / 2 / (1 + NU), E * NU / (1 - NU * NU)
    sig = lambda w: 2 * mu * Sym(Grad(w)) + lam * Trace(Sym(Grad(w))) * Id(2)
    a = BilinearForm(InnerProduct(sig(u), Sym(Grad(v))) * dx).Assemble()
    l = LinearForm(v[1] * ds("load")).Assemble()
    g = GridFunction(fes)
    g.vec.data = a.mat.Inverse(fes.FreeDofs(), inverse="sparsecholesky") * l.vec
    xs = np.linspace(1e-6, xr - 1e-6, 3001)
    sy = np.array([sig(g)[1, 1](mesh(x, 1e-9)) for x in xs])
    return dict(order=order, Kt_net=float(sy[-1] * (M.W - 2 * M.C) / M.W), ligament_force=float(np.trapezoid(sy, xs)),
                applied=M.W / 2)


def neuber_deep(a, rho):
    q = np.sqrt(a / rho)
    return 2 * (a / rho + 1) * q / ((a / rho + 1) * np.arctan(q) + q)


def main():
    rows = [run_level(h) for h in LEVELS]
    rows[0]["anchors"] = dict(plane_stress_2d=[kt_2d(p) for p in (2, 3, 4)],
                              neuber_deep=float(neuber_deep((M.W - 2 * M.C) / 2, R_NOTCH)))
    json.dump(rows, open(os.path.join(OUT, "ccx_ngsolve.json"), "w"), indent=1)
    return rows


def report(rows=None):
    rows = rows or json.load(open(os.path.join(OUT, "ccx_ngsolve.json")))
    L = ["# Toolchain check: NGSolve P2 vs CalculiX C3D10 -> critical-plane search", "",
         "Ezeh & Susmel (2020) intermediate notch (U, r = 1 mm), 3D 1/8 model, PLA E = 3476 MPa, nu = 0.33. Identical "
         "straight-sided tetrahedra and mid-edge nodes for both solvers (same discrete space, exact stiffness "
         "integration), nodal stresses averaged over elements in both. Net-section amplitude 10 MPa. 4 CPU cores.", "",
         "## 1. FE agreement and mesh convergence", "",
         "| notch h [mm] | vertices | DOF NGSolve / ccx | NGSolve solve [s] | ccx solve [s] | ccx write+read [s] | "
         "max dU / U | max dS / S | Kt_net surface NGS / ccx | Kt_net mid-plane NGS / ccx |", "|" + "---|" * 10]
    for r in rows:
        L.append("| %.2f | %d | %d / %d | %.2f | %.2f | %.2f | %.1e | %.1e | %.3f / %.3f | %.3f / %.3f |" % (
            r["h"], r["nv"], r["ndof_ngs"], r["ndof_ccx"], r["t_ngs"], r["t_ccx"], r["t_ccx_io"], r["dU_max_rel"],
            r["dS_max_rel"], r["Kt_surf_ngs"], r["Kt_surf_ccx"], r["Kt_mid_ngs"], r["Kt_mid_ccx"]))
    an = rows[0].get("anchors")
    if an:
        L += ["", "Anchors: 2D plane stress on the curved geometry (notch edge h = 0.03 mm), Kt_net = " + ", ".join(
            "%.4f (p = %d)" % (q["Kt_net"], q["order"]) for q in an["plane_stress_2d"]) +
              "; ligament force / applied = %.5f. Neuber deep hyperbolic notch (a = 3, rho = 1 mm): %.3f. "
              "3D brackets plane stress: mid-plane (constrained) above, free surface below." % (
                  an["plane_stress_2d"][-1]["ligament_force"] / an["plane_stress_2d"][-1]["applied"], an["neuber_deep"]),
              "", "**Note:** the earlier elasticity_rhmp EFFM study (fakasht10/codex) reported Kt_net 4.18 / 4.02 for this notch; an "
              "independent audit traced that to its root-stress recovery (nodal force / line weight). The anchors here hold."]
    L += ["", "frd stress components (checked against cplane order): " + ", ".join(rows[-1]["frd_components"][4:10]), "",
          "## 2. Critical-plane results on the free-surface nodes (finest mesh)", "",
          "| R | criterion | max damage param. NGS [MPa] | ccx [MPa] | rel. diff | max field diff / max | same hot-spot node | "
          "normal diff [deg] | hot-spot normal (x, y, z) | tau_a / sn_a / sn_m [MPa] | search [s] | screened [s] (points) |",
          "|" + "---|" * 12]
    r = rows[-1]
    for k, c in r["cplane"].items():
        R, cname = k.split("|")
        L.append("| %s | %s | %.4f | %.4f | %.1e | %.1e | %s | %.2f | (%.2f, %.2f, %.2f) | %.2f / %.2f / %.2f | %.3f | %.3f (%d of %d)%s |" % (
            R[2:], cname, c["vmax_ngs"], c["vmax_ccx"], c["dv"], c["field_drel"], "yes" if c["same_node"] else "no",
            c["dnormal_deg"], *c["normal"], c["tau_a"], c["sn_a"], c["sn_m"], c["t"], c["t_screen"], c["n_screen"],
            r["n_surf"], "" if c["screen_ok"] else " MISSED"))
    L += ["", "Mesh convergence of the hot-spot damage parameter (NGSolve, Findley k = 0.3, R = -1): " + ", ".join(
        "h %.2f: %.4f" % (q["h"], q["cplane"]["R=-1|Findley k=0.3"]["vmax_ngs"]) for q in rows) + " MPa.", "",
          "![fields](ccx_ngsolve.png)"]
    open(os.path.join(OUT, "CCX_NGSOLVE.md"), "w").write("\n".join(L) + "\n")
    import matplotlib.tri as mtri
    fig, ax = plt.subplots(1, 3, figsize=(13.5, 4.3), gridspec_kw=dict(width_ratios=[1.1, 1, 1]), layout="constrained")
    h = [q["h"] for q in rows]
    ax[0].semilogx(h, [q["Kt_surf_ngs"] for q in rows], "o-", label="3D free surface, NGSolve P2")
    ax[0].semilogx(h, [q["Kt_surf_ccx"] for q in rows], "x--", label="3D free surface, ccx C3D10")
    ax[0].semilogx(h, [q["Kt_mid_ngs"] for q in rows], "s-", label="3D mid-plane, NGSolve P2")
    ax[0].semilogx(h, [q["Kt_mid_ccx"] for q in rows], "+--", label="3D mid-plane, ccx C3D10")
    if an:
        ax[0].axhline(an["plane_stress_2d"][-1]["Kt_net"], color="0.3", lw=1, ls=":", label="2D plane stress (p = 4)")
        ax[0].axhline(an["neuber_deep"], color="0.6", lw=1, ls="-.", label="Neuber deep notch")
    ax[0].set(xlabel="notch mesh size h [mm]", ylabel=r"$K_{t,net}$", title="(a) notch-root Kt: both solvers, anchors")
    ax[0].invert_xaxis()
    ax[0].legend(fontsize=6.5, loc="upper center", bbox_to_anchor=(0.5, -0.2), ncol=2)
    ax[0].grid(alpha=0.3, which="both")
    pl = rows[-1]["_plot"]
    X, F, tri, xr = np.array(pl["X"]), np.array(pl["F"]), np.array(pl["tri"]), pl["xr"]
    top = tri[np.all(np.abs(X[tri, 2] - M.T / 2) < 1e-9, axis=1)]
    vmax = np.nanmax(F)
    tp = mtri.Triangulation(X[:, 0], X[:, 1], top)
    ax[1].tripcolor(tp, np.nan_to_num(F), shading="gouraud", cmap="viridis", vmin=0, vmax=vmax)
    ax[1].set(xlim=(0, 6), ylim=(0, 3.5), aspect="equal", xlabel="x [mm]", ylabel="y [mm]",
              title="(b) free face z = t/2 (zoom)")
    cx = xr + R_NOTCH
    notch = tri[np.all(np.hypot(X[tri, 0] - cx, X[tri, 1]) < R_NOTCH + 1e-6, axis=1)]
    th = np.degrees(np.arctan2(X[:, 1], -(X[:, 0] - cx)))                # 0 deg = root, 90 deg = notch flank start
    tn = mtri.Triangulation(th, X[:, 2], notch)
    im = ax[2].tripcolor(tn, np.nan_to_num(F), shading="gouraud", cmap="viridis", vmin=0, vmax=vmax)
    i = np.nanargmax(F)
    ax[2].plot(th[i], X[i, 2], "r*", ms=14, clip_on=False, zorder=5, label="hot spot (root, mid-plane)")
    ax[2].set(xlabel="angle from notch root [deg]", ylabel="z [mm] (0 = mid-plane)", xlim=(0, 90), ylim=(0, M.T / 2),
              title="(c) notch surface, unrolled")
    ax[2].legend(fontsize=7, loc="upper right")
    fig.colorbar(im, ax=ax[1:], label="SWT-stress parameter, R = 0 [MPa]", shrink=0.85)
    _save(fig, OUT, "ccx_ngsolve")


if __name__ == "__main__":
    report(main())
