# ECP step 2: specimens rebuilt in NGSolve, ECP factors from `cplane`, compared with Chiocca & Frendo

**Summary.**
* **The independent chain reproduces their factors.** At r_c = 0.2 mm, 11 of the 13 cylindrical series agree within 2.0 % and the flat plates within 7.3 % (ours lower). Kt anchors agree within 2-5 %.
* **Two 90 deg R = -1 series differ, and ours is the global maximum.** V-notch r = 0.2 is 1.26x theirs and r = 4 is 1.09x. Our plane is confirmed by a 200k-normal brute force. Their value lies below even the classic max-shear-plane FS for r = 0.2, and equals it for r = 4, which is consistent with a local optimum in their search.
* **Predictions survive the swap of FE code.** With their calibration procedure on our factors, the prevision E_N is +0.07 ± 0.42 with 74 % within a factor of 3 (theirs: -0.00 ± 0.44, 73 %). The calibrated r_c is the same, 0.220 mm. So the ECP result is not an artefact of their FE models.
* **The weak point is non-proportional loading, not FE.** With the plain MCC shear amplitude, the 90 deg series are non-conservative, with mean E_N -0.51 (ours) and -0.68 (theirs). With their MPHP factors they become conservative, +0.41 and +0.25. Neither amplitude measure centres the 90 deg data. Our higher global-max factors also push the 90 deg series further conservative.

**Checklist.**
* **Theory.** FS = (τ_a / G)(1 + k σ_n,max / σ_y) on averaged tensors, with k = 0.4 and σ_y = 606.2. The extended plane maximises the whole parameter. Uniaxial rows use the max-shear plane at unit load, as their rows do. The MCC shear amplitude is exact (Welzl).
* **Physics.** Linear elastic En3B (E = 200 GPa, ν = 0.3).
  * Bars: 2D axisymmetric (u_r, u_z) plus a separate torsion problem (u_θ), on net nominal stress.
  * Plates: 3D, clamped bottom, uniform traction on top (free rotation, so the single-edge notches bend), mid-plane symmetry, F = 1 N on a 150 mm² gross area.
  * Hot spot: on the root line, maximum τ_a. The V-notch peak is at mid-plane, 1.38× the surface value.
* **Numerics.** P2 elements on straight-sided simplices.
  * Local h: bars 0.02 → 0.01 mm (root ≤ r_n / 15); plates 0.08 → 0.05 mm; V-notch plate 0.03 → 0.02 mm in a 0.35 mm ball, with r_n / 6 at the root.
  * Solvers: direct sparse Cholesky, or CG + BDDC above 250k DOF.
  * Averages: quadrature-point integral means, with the 3D sphere mirrored across the mid-plane. Below the quadrature spacing the point value is used.
  * Mesh convergence: averaged tensors change by ≤ 2.5 % between the two levels (table below).
  * Open: Kt sits 2–5 % above their FE values (the BC or nominal definition may differ). The flat V-notch curve dips 12–17 % below their 6th-order fit at r_c = 0.05–0.15 mm, but is mesh-converged at the calibrated 0.2 mm.
* **Visualisation.** `ecp_fe_factor_curves` compares the factor curves for 8 series. `ecp_fe_end_to_end` shows the bias and scatter per loading family, for both factor sets and both amplitude measures.

## Mesh / Kt anchors
| specimen | Kt paper (FE) | Kt here, coarse -> fine | DOF fine | max change of averaged tensor coarse->fine |
|---|---|---|---|---|
| V_notch_02 | 3.86 / 2.02 | 3.786 / 2.035 -> 3.786 / 2.036 | 53728 | 0.26 % |
| V_notch_12 | 1.776 / 1.317 | 1.777 / 1.269 -> 1.777 / 1.269 | 40975 | 0.15 % |
| V_notch_4 | 1.278 / 1.121 | 1.273 / 1.096 -> 1.273 / 1.096 | 37384 | 0.03 % |
| V_notch | 16.52 | 17.335 -> 16.916 | 264453 | 0.69 % |
| U_notch | 6.21 | 6.413 -> 6.406 | 193636 | 1.26 % |
| Hole_8 | 3.45 | 3.542 -> 3.558 | 106690 | 0.75 % |
| Hole_3_5 | 3.13 | 3.167 -> 3.184 | 196420 | 2.49 % |

Bars: tension / torsion Kt on the net section. Plates: gross Kt (max principal stress at the hot spot x 150 mm^2).

## Our ECP factors / theirs (same convention, MCC amplitude, no corrector)
| series | r_c = 0 | 0.1 | 0.2 | 0.3 mm |
|---|---|---|---|---|
| Hole 8 mm, R=-1 | 0.986 | 0.967 | 0.974 | 0.978 |
| V-notch (flat), R=-1 | 0.962 | 0.811 | 0.937 | 1.020 |
| Hole 3.5 mm, R=-1 | 0.986 | 0.955 | 0.954 | 0.962 |
| Hole 3.5 mm, R=0.1 | 0.985 | 0.954 | 0.953 | 0.961 |
| Hole 8 mm, R=0.1 | 0.984 | 0.972 | 0.973 | 0.974 |
| U-notch (flat), R=-1 | 1.003 | 0.978 | 0.974 | 0.971 |
| U-notch (flat), R=0.1 | 1.003 | 0.977 | 0.974 | 0.970 |
| V-notch (flat), R=0.1 | 0.962 | 0.801 | 0.927 | 1.018 |
| V-notch r=0.2, 0 deg, R=-1 | 1.012 | 0.974 | 1.009 | 0.997 |
| V-notch r=0.2, 0 deg, R=0 | 1.028 | 1.021 | 1.020 | 1.012 |
| V-notch r=0.2, 90 deg, R=-1 | 1.011 | 1.008 | 1.264 | 1.373 |
| V-notch r=0.2, 90 deg, R=0 | 0.991 | 0.993 | 0.999 | 0.997 |
| V-notch r=1.25, 0 deg, R=-1 | 1.020 | 1.000 | 0.995 | 0.994 |
| V-notch r=1.25, 0 deg, R=0 | 1.032 | 1.020 | 1.013 | 1.002 |
| V-notch r=1.25, 90 deg, R=-1 | 1.022 | 0.997 | 0.994 | 1.012 |
| V-notch r=1.25, 90 deg, R=0 | 1.010 | 0.996 | 0.995 | 0.993 |
| V-notch r=4, 0 deg, R=-1 | 1.009 | 1.006 | 1.006 | 1.005 |
| V-notch r=4, 0 deg, R=0 | 1.017 | 1.015 | 1.014 | 1.012 |
| V-notch r=4, 90 deg, R=-1 | 1.101 | 1.094 | 1.090 | 1.085 |
| V-notch r=4, 90 deg, R=0 | 1.005 | 1.003 | 1.003 | 1.002 |
| V-notch r=0.2, torsion, R=-1 | 0.997 | 0.991 | 0.996 | 0.994 |

## End to end: calibrate on V-notch + hole 8 (R = -1), predict the 19 other series
Cells: mean E_N / std E_N / share within factor 3. r_c = 0.2 mm except the 'r_c opt' column.

| factors | r_c opt [mm] | all 19 | all 19 at r_c opt | flat axial | in-phase | 90 deg | torsion |
|---|---|---|---|---|---|---|---|
| their factors (ANSYS), MPHP corr. | 0.220 | -0.00 / 0.44 / 73 % | +0.01 / 0.44 / 74 % | -0.03 / 0.31 / 91 % | -0.25 / 0.40 / 62 % | +0.25 / 0.49 / 61 % | +0.11 / 0.40 / 67 % |
| our factors (NGSolve), MPHP corr. | 0.220 | +0.07 / 0.42 / 74 % | +0.07 / 0.43 / 74 % | -0.04 / 0.31 / 91 % | -0.15 / 0.39 / 72 % | +0.41 / 0.38 / 56 % | +0.18 / 0.40 / 56 % |
| our factors (NGSolve), MCC only | 0.220 | -0.19 / 0.43 / 70 % | -0.18 / 0.44 / 73 % | -0.04 / 0.31 / 91 % | -0.15 / 0.39 / 72 % | -0.51 / 0.44 / 41 % | +0.18 / 0.40 / 56 % |
| their factors (ANSYS), MCC only | 0.220 | -0.26 / 0.50 / 66 % | -0.24 / 0.50 / 67 % | -0.03 / 0.31 / 91 % | -0.25 / 0.40 / 62 % | -0.68 / 0.55 / 34 % | +0.11 / 0.40 / 67 % |

Wall time of the whole study: 537 s.

