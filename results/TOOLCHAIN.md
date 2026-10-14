# Fatigue toolchain: FE stress -> fast critical-plane search -> durability

**Recommendation.** Use NGSolve as the primary solver: it runs in-process, uses our complexes, makes superposition
and anisotropic C trivial, and was 1.5–2.4× faster than ccx here. Keep CalculiX as an independent cross-check, and
for decks that need contact or nonlinear materials. Both feed `cplane.py` through the same (points × 6) stress
arrays. Do not adopt AlejandroBSantana/Abaqus-Fatigue-Fracture-Tools; reasons are below.

## Evidence (this folder)

| check | result | report |
|---|---|---|
| kernel vs closed forms | Findley uniaxial and torsion exact to 1e-8; MCC exact on circle and collinear paths, ellipse within 1e-3 (discrete sampling) | `tests/test_cplane.py` (14 tests) |
| kernel vs brute force | prop ≤ 6e-7 max rel. error (reference: kernel at 8000 normals; independent audit: isotropic ≤ 1.2e-9, anisotropic with a crease-circle search ≤ 1e-9 below a 300k-normal + crease reference); path ≤ 1.5e-4 vs a 200 000-normal grid without polish (part of it is the reference's own error) | `CPLANE_SPEED.md` |
| speed | 4.9e4 points/s full field (1e6 points in 20 s); 1e6 points with hot-spot screening in 2.2 s (synthetic log-normal field: only 4096 points needed an exact search; FE fields differ); naive loop 15 points/s; dense 1° grid ~55× slower and ~300× less accurate | `CPLANE_SPEED.md` |
| NGSolve P2 vs ccx C3D10, same mesh | identical DOF; dU 2.8e-6, dS 1.4e-5 (.frd has 5 significant digits); same hot-spot node and normal; damage parameter within 1.2e-5 | `CCX_NGSOLVE.md` |
| Kt anchors, U r = 1 notch | 3D mid-plane 2.606 / surface 2.212 (converged); 2D plane stress 2.475; Neuber 2.34 | `CCX_NGSOLVE.md` |
| solve time, 86k DOF | NGSolve 2.9–5.4 s (sparse Cholesky) vs ccx 6.9–10 s (SPOOLES) + 0.3 s file I/O | `CCX_NGSOLVE.md` |

**Open issue.** `EFFM.md` reports Kt_net = 4.18 (2D) and 4.02 (3D) for the same notch. That is 1.7× the anchors
above, so its Kt column or normalisation should be checked.

## CalculiX: viable, with these gotchas (all handled in `ccx_io.py`)
* Input fields are cut at 20 characters, silently or with misleading errors (for example "Poisson > 0.5"). Write
  numbers as `%.12e`.
* .frd results carry 5 significant digits. This is enough for fatigue; use `*EL PRINT` (.dat) for more.
* .frd STRESS order is SXX SYY SZZ SXY SYZ SZX, which equals the cplane order. Abaqus and .dat order is
  11 22 33 12 13 23; convert with `cplane.from_abaqus`.
* C3D10 needs positive orientation and the ccx mid-edge order 5:(1,2) 6:(2,3) 7:(3,1) 8:(1,4) 9:(2,4) 10:(3,4).
  Pressure faces are F1..F4.
* **Run the ccx equation solver single-threaded** (`CCX_NPROC_EQUATION_SOLVER=1`): multi-threaded SPOOLES in
  ccx 2.21 intermittently returns a wrong solution (1 of 12 identical runs, 35 % off). This, not vertex
  numbering, explains the one-off NGSolve/ccx mismatch seen at the h = 0.2 mesh.
* Netgen and NGSolve can renumber vertices: always assert coordinates before mapping results between codes. This
  check caught a real mismatch here.

## Abaqus-Fatigue-Fracture-Tools: reviewed, not usable for this work
* **It is not a critical-plane method.** `Criterion_Dang_Van.py` evaluates Dang Van (mesoscopic Tresca plus
  hydrostatic stress) from principal stresses; there is no plane search. The other scripts are Goodman on von Mises
  or on max principal stress.
* **Physics is wrong for printed polymers.** Dang Van rests on elastic shakedown of metal grains. That does not apply
  to amorphous PLA or ABS, nor to raster or inter-layer failure.
* **Numerics:**
  * Badoiu–Clarkson enclosing ball with a fixed 1000 iterations from the mean. This is only ε-approximate (worst-case
    bound about 3 % on the radius for 1000 steps; typical errors are smaller), and the residual-stress centre error feeds straight into the index.
  * Pure-Python loops over integration points × increments.
  * Increments above 900 are skipped silently.
* **Platform:** Abaqus 6.14 `odbAccess`, Python 2.7. ODB, instance, set and material (Ti-6Al-4V) are all hard-coded.
  Output is the element-centroid worst value only.
* **No licence file,** so the code cannot be reused.
* **Worth keeping:** the idea of a minimum enclosing ball on the deviatoric path. If metal AM needs Dang Van or
  Papadopoulos, add it to `cplane` with an exact enclosing-ball method (Welzl), not B-C.

## Theory / physics / numerics / visualisation checklist (as built)
* **Theory.** Stress-based critical plane (Findley; MWCM with ρ_eff; elastic SWT for mode-I inter-layer planes).
  τ_a is the minimum circumscribed circle radius on each plane. Anisotropy enters as strength(|n·d|), so it changes
  which plane is critical.
* **Physics.** Linear-elastic, constant amplitude, proportional or one general cycle. Not yet covered: variable
  amplitude (rainflow on each plane), self-heating, viscoelastic effects, calibration of the criterion constants
  (k, τ_A/σ_A, strength polynomial: these need torsion or multiaxial PLA/ABS data).
* **Numerics.**
  1. Per-point Numba kernel.
  2. Coarse Fibonacci hemisphere: 800 normals for proportional loading, 400 for paths.
  3. Refine the 4 best coarse local maxima (8-neighbour graph). Findley has near-tied peaks tens of degrees apart,
     so a top-q-values rule fails.
  4. Stencil Newton with a golden-angle-rotated compass fallback.
  5. For paths, a polish stage on MCC creases.
  6. Upper-bound hot-spot screening.
* **Visualisation.** Time/accuracy log-log plots; Kt convergence with anchors; damage maps on the free face and on
  the unrolled notch surface.
