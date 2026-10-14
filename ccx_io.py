"""Minimal CalculiX (ccx) bridge: write a C3D10 linear-elastic deck from a tetrahedral mesh, run ccx, read the .frd.

    deck = write_c3d10_inp(path, X, T, E, nu, fix={1: nodes_x0, 2: ..., 3: ...}, pressure_faces=[(el, face), ...])
    run_ccx(jobdir, "job")                    # -> job.frd
    nodes, U, S = read_frd("job.frd")         # S in [xx, yy, zz, xy, yz, zx] (cplane order)

Mid-edge nodes are created on straight edges, so a C3D10 model is the same discrete space as an NGSolve VectorH1
order-2 field on the identical straight-sided mesh (and the 4-point rule of C3D10 integrates its stiffness exactly).
ccx reads at most 20 characters per field: numbers are written as %.12e.
ccx node / face conventions (ccx manual, C3D10): corners 1-4, mid-edges 5:(1,2) 6:(2,3) 7:(3,1) 8:(1,4) 9:(2,4) 10:(3,4);
faces F1:(1,2,3) F2:(1,4,2) F3:(2,4,3) F4:(3,4,1).  .frd values are written as E12.5 (5 significant digits).
"""
import os
import subprocess
import time

import numpy as np

C3D10_EDGES = ((0, 1), (1, 2), (2, 0), (0, 3), (1, 3), (2, 3))
C3D10_FACES = ((0, 1, 2), (0, 3, 1), (1, 3, 2), (2, 3, 0))


def orient(X, T):
    """Positive-volume corner ordering (ccx stops on negative Jacobians)."""
    T = T.copy()
    a, b, c, d = (X[T[:, i]] for i in range(4))
    neg = np.einsum("ij,ij->i", np.cross(b - a, c - a), d - a) < 0
    T[neg, 1], T[neg, 2] = T[neg, 2].copy(), T[neg, 1].copy()
    return T


def quadratic(X, T):
    """Add straight mid-edge nodes. Returns X10 (n10, 3), T10 (ne, 10) in ccx order; corner ids unchanged."""
    E = np.sort(np.concatenate([T[:, [i, j]] for i, j in C3D10_EDGES]), axis=1)
    uniq, inv = np.unique(E, axis=0, return_inverse=True)
    mid = len(X) + inv.reshape(len(C3D10_EDGES), len(T)).T
    X10 = np.vstack([X, 0.5 * (X[uniq[:, 0]] + X[uniq[:, 1]])])
    return X10, np.hstack([T, mid])


def write_c3d10_inp(path, X10, T10, E, nu, fix, pressure_faces, pressure=-1.0, threads=4):
    """fix: {dof: node ids (0-based)}; pressure_faces: (element 0-based, face 0..3); ccx pressure > 0 pushes in."""
    with open(path, "w") as f:
        f.write("*HEADING\ncplane bridge: C3D10 linear elastic\n*NODE, NSET=NALL\n")
        np.savetxt(f, np.column_stack([np.arange(1, len(X10) + 1), X10]), fmt="%d, %.12e, %.12e, %.12e")
        f.write("*ELEMENT, TYPE=C3D10, ELSET=EALL\n")
        np.savetxt(f, np.column_stack([np.arange(1, len(T10) + 1), T10 + 1]), fmt="%d", delimiter=", ")
        for dof, nodes in fix.items():
            f.write(f"*NSET, NSET=FIX{dof}\n")
            np.savetxt(f, np.asarray(nodes) + 1, fmt="%d")
        f.write("*MATERIAL, NAME=M\n*ELASTIC\n%.12e, %.12e\n*SOLID SECTION, ELSET=EALL, MATERIAL=M\n" % (E, nu))
        f.write("*STEP\n*STATIC, SOLVER=SPOOLES\n*BOUNDARY\n")
        for dof in fix:
            f.write(f"FIX{dof}, {dof}, {dof}\n")
        f.write("*DLOAD\n")
        for e, k in pressure_faces:
            f.write("%d, P%d, %.12e\n" % (e + 1, k + 1, pressure))
        f.write("*NODE FILE\nU\n*EL FILE\nS\n*END STEP\n")
    return path


def run_ccx(jobdir, job, threads=4):
    # The equation solver runs single-threaded: multi-threaded SPOOLES in ccx 2.21 returned a wrong solution in
    # 1 of 12 identical runs (35 % off) on a 24k-DOF C3D10 model; 1 thread: 12 / 12 identical. Stiffness assembly may
    # still use `threads`.
    env = dict(os.environ, OMP_NUM_THREADS=str(threads), CCX_NPROC_EQUATION_SOLVER="1",
               CCX_NPROC_STIFFNESS=str(threads))
    t = time.perf_counter()
    r = subprocess.run(["ccx", "-i", job], cwd=jobdir, env=env, capture_output=True, text=True)
    dt = time.perf_counter() - t
    if r.returncode != 0 or "ERROR" in r.stdout:
        raise RuntimeError(r.stdout[-3000:] + r.stderr[-2000:])
    return dt, r.stdout


def read_frd(path):
    """Nodal DISP and STRESS (extrapolated from integration points and averaged over elements by ccx)."""
    blocks, cur, name = {}, None, None
    with open(path) as f:
        for line in f:
            if line.startswith(" -4"):
                name = line.split()[1]
                cur = blocks.setdefault(name, {})
            elif line.startswith(" -1") and name is not None:
                node = int(line[3:13])
                vals = [float(line[13 + 12 * i:25 + 12 * i]) for i in range((len(line.rstrip("\n")) - 13) // 12)]
                cur[node] = vals
            elif line.startswith(" -3"):
                name = None
    U, S = blocks["DISP"], blocks["STRESS"]
    nodes = np.array(sorted(S))
    Ua = np.array([U[n][:3] for n in nodes])
    Sa = np.array([S[n][:6] for n in nodes])        # frd STRESS: SXX SYY SZZ SXY SYZ SZX == cplane order
    return nodes - 1, Ua, Sa
