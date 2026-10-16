# ECP (Chiocca & Frendo 2024) reimplemented and stress-tested on its own data

**Summary.**
* **Reproduction is exact.** The Python port reproduces the repository's MATLAB outputs to 4e-15 relative: calibration T = 1.2025, prevision T = 1.7666 and all 19 per-series T values. The calibration residual curve matches to 2e-13.
* **The headline prevision T = 1.77 is a self-fit, but here it happens to be honest.** Like the paper's Table 4, it refits a new power law (m = -0.185, against the calibrated m = -0.230). Scored against the calibrated curve, as a prediction must be, the same T formula gives 1.80.
* **The near-zero multiaxial bias is a cancellation.** In-phase series are non-conservative (mean E_N -0.25) and 90 deg series conservative (+0.25). The 90 deg result rests on the MPHP corrector: with the plain MCC shear amplitude (what `cplane` and most codes use) the 90 deg mean is -0.68, about a factor 5 non-conservative.
* **Predicted lives are still good for a one-parameter model.** On the 19 prevision series, E_N has mean -0.00 and std 0.44; 73 % fall within a factor of 3 (the paper reports 74 % on all data).
* **Averaging is what makes it work.** At r_c = 0 (hot spot, no averaging) the rms E_N is 4.09 decades, against 0.43 at r_c = 0.2 mm. The oracle r_c, chosen on the prevision data, is 0.196 mm with rms 0.43. So the two-series calibration picks a near-optimal radius.
* **r_c is robust to the calibration choice.** Re-running the paper's 26 combinations (Table 3) gives r_c = 0.199 ± 0.017 mm (paper: 0.204 ± 0.016); 20 of 26 radii match the paper within 0.01 mm.

## A. Reproduction of the repository outputs (r_c = 0.2 mm grid point)
| quantity | Python | repository (MATLAB) |
|---|---|---|
| calibration T | 1.20248 | 1.20248 |
| calibration c | 0.0259125 | 0.0259125 |
| calibration m | -0.229523 | -0.229523 |
| prevision T (self-fit) | 1.76658 | 1.76658 |
| prevision m (self-fit) | -0.185353 | -0.185353 |

| prevision series | T (self-fit), Python | repository |
|---|---|---|
| Axial_Hole_3_5_R_-1 | 1.06018 | 1.06018 |
| Axial_Hole_3_5_R_0.1 | 1.15048 | 1.15048 |
| Axial_Hole_8_R_0.1 | 1.10605 | 1.10605 |
| Axial_U_notch_R_-1 | 1.07425 | 1.07425 |
| Axial_U_notch_R_0.1 | 1.09810 | 1.09810 |
| Axial_V_notch_R_0.1 | 1.39196 | 1.39196 |
| Multi_V_notch_02_Ph_0_R_-1 | 1.27086 | 1.27086 |
| Multi_V_notch_02_Ph_0_R_0 | 1.13169 | 1.13169 |
| Multi_V_notch_02_Ph_90_R_-1 | 1.32555 | 1.32555 |
| Multi_V_notch_02_Ph_90_R_0 | 1.26201 | 1.26201 |
| Multi_V_notch_12_Ph_0_R_-1 | 1.10065 | 1.10065 |
| Multi_V_notch_12_Ph_0_R_0 | 1.06643 | 1.06643 |
| Multi_V_notch_12_Ph_90_R_-1 | 1.18887 | 1.18887 |
| Multi_V_notch_12_Ph_90_R_0 | 1.20242 | 1.20242 |
| Multi_V_notch_4_Ph_0_R_-1 | 1.20092 | 1.20092 |
| Multi_V_notch_4_Ph_0_R_0 | 1.25179 | 1.25179 |
| Multi_V_notch_4_Ph_90_R_-1 | 1.16524 | 1.16524 |
| Multi_V_notch_4_Ph_90_R_0 | 1.09504 | 1.09504 |
| Torsion_V_notch_02_R_-1 | 1.19439 | 1.19439 |

Optimal r_c on the repository grid (0-0.5 mm, 500 points): 0.2194 mm. MAIN.m evaluates at the fixed value 0.2.

## B. Predictive error against the calibrated curve, r_c = 0.2 mm
E_N = log10(N_f / N_exp), paper Eq. 9; > 0 means conservative. All rows are predictions except where a group includes the two calibration series (marked 'all').

| group | n | mean E_N | std E_N | % in x2 | % in x3 | % in x5 |
|---|---|---|---|---|---|---|
| prevision (19 series) | 147 | -0.003 | 0.436 | 45 | 73 | 89 |
| flat axial, prevision | 57 | -0.027 | 0.308 | 63 | 91 | 98 |
| flat axial, all (paper Fig. 11a) | 69 | -0.022 | 0.285 | 70 | 93 | 99 |
| cylindrical multiaxial + torsion (paper Fig. 11b) | 90 | +0.013 | 0.502 | 33 | 62 | 83 |
|   in-phase (0 deg) | 40 | -0.250 | 0.400 | 35 | 62 | 88 |
|   out-of-phase (90 deg) | 41 | +0.249 | 0.494 | 29 | 61 | 76 |
|   torsion | 9 | +0.107 | 0.397 | 44 | 67 | 100 |

Paper Fig. 11 (FS ECP): uniaxial mean -0.022, std 0.285; multiaxial mean +0.013, std 0.492.

**Scatter on the paper's T scale** (T = exp(2 · 1.282 · s · |m|), s = std of ln N):

| definition | T |
|---|---|
| calibration data vs calibration curve (paper: 1.20) | 1.202 |
| prevision vs its own refitted curve (paper Table 4: 1.77) | 1.767 |
| prevision vs calibration curve, bias removed | 1.806 |
| prevision vs calibration curve, rms (bias kept) | 1.802 |

Per-series bias (mean ± std of E_N) is in `ecp_series_bias.png`. Every series has a T near 1.1-1.4 about its own line (Table 4), but the series lines are offset from each other. The offsets are what a designer would actually incur.

## C. Sensitivities (19 prevision series)
| variant | mean E_N | std | % x2 | % x3 | % x5 |
|---|---|---|---|---|---|
| r_c = 0 (hot spot, no averaging) | -1.139 | 3.938 | 8 | 10 | 14 |
| r_c = 0.2 mm (paper) | -0.003 | 0.436 | 45 | 73 | 89 |
| r_c = 0.219 mm (calibration optimum) | +0.011 | 0.443 | 47 | 74 | 88 |
| r_c = 0.196 mm (oracle on prevision data; not a prediction) | -0.006 | 0.436 | 48 | 73 | 90 |
| r_c = 0.2, regression ln N on ln FS (equivalent m = -0.244 vs -0.230) | -0.026 | 0.422 | 48 | 74 | 90 |

90 deg out-of-phase series only (the corrector multiplies FS by the MPHP / MCC ratio, 1.46-1.77):

| shear amplitude | mean E_N | std | % x2 | % x3 | % x5 |
|---|---|---|---|---|---|
| MPHP (paper) | +0.249 | 0.494 | 29 | 61 | 76 |
| MCC (no corrector) | -0.684 | 0.555 | 22 | 34 | 56 |

## D. Paper Table 3 re-run (N_f > 3e4, r_c grid 0-0.7 mm in 0.001 steps)
| combination (geometry, R) | r_c paper | r_c here | 1/|m| paper | 1/|m| here |
|---|---|---|---|---|
| U-1 V-1 | 0.19 | 0.186 | 3.62 | 3.77 |
| U.1 V.1 | 0.18 | 0.171 | 3.12 | 3.31 |
| H8-1 V-1 | 0.20 | 0.219 | 4.29 | 4.17 |
| H8.1 V.1 | 0.23 | 0.208 | 3.39 | 3.40 |
| H3.5.1 U.1 | 0.19 | 0.190 | 3.85 | 4.05 |
| H3.5-1 U-1 | 0.19 | 0.195 | 4.39 | 4.53 |
| H3.5-1 V-1 | 0.19 | 0.192 | 4.48 | 4.31 |
| H3.5.1 V.1 | 0.19 | 0.195 | 3.22 | 3.56 |
| H3.5-1 V-1 H3.5.1 V.1 | 0.20 | 0.189 | 3.69 | 3.92 |
| H3.5-1 V-1 H3.5.1 U.1 | 0.21 | 0.212 | 4.18 | 4.20 |
| H3.5-1 V-1 H3.5.1 V.1 | 0.20 | 0.189 | 3.69 | 3.92 |
| U-1 V-1 H8-1 H3.5-1 | 0.20 | 0.201 | 4.10 | 4.21 |
| U.1 V.1 H8.1 H3.5.1 | 0.20 | 0.192 | 3.52 | 3.53 |
| H8-1 V-1 H8.1 V.1 | 0.21 | 0.206 | 3.73 | 3.80 |
| H3.5-1 V-1 H8.1 V.1 | 0.21 | 0.195 | 3.80 | 3.83 |
| H3.5-1 V-1 H8.1 U.1 | 0.23 | 0.227 | 4.42 | 3.95 |
| U-1 V-1 U.1 V.1 | 0.18 | 0.172 | 3.40 | 3.60 |
| U-1 V-1 U.1 H8.1 | 0.24 | 0.231 | 3.68 | 3.47 |
| U-1 V-1 U.1 H3.5.1 | 0.21 | 0.214 | 3.62 | 3.75 |
| H8-1 V-1 U.1 V.1 | 0.19 | 0.187 | 3.45 | 3.71 |
| H8-1 V-1 U.1 H8.1 | 0.23 | 0.234 | 4.17 | 3.85 |
| H8-1 V-1 U.1 H3.5.1 | 0.22 | 0.221 | 3.95 | 4.09 |
| H3.5-1 V-1 U.1 V.1 | 0.19 | 0.174 | 3.57 | 3.78 |
| U-1 V-1 V.1 H8.1 | 0.20 | 0.194 | 3.50 | 3.53 |
| U-1 V-1 V.1 H3.5.1 | 0.19 | 0.187 | 3.39 | 3.66 |
| H8-1 V-1 V.1 H3.5.1 | 0.20 | 0.199 | 3.57 | 3.89 |

Max |Δr_c| = 0.022 mm, max |Δ(1/m)| = 0.47. Table rows 9 and 11 are identical in the paper; both are kept as printed.

## Theory / physics / numerics / visualisation checklist
* **Theory.** Fatemi-Socie on the averaged tensors: FS = (Δγ/2)(1 + k σ_n,max / σ_y), k = 0.4, σ_y = 606.2 MPa. Linear elasticity, so FS(F) = (F Δγ₁ / 2)(1 + k F σ_n1 / σ_y) from unit-load fields.
* **Physics.** En3B is a ductile low-carbon steel with shear-dominated crack initiation, which is why FS suits it (the paper's SWT variant scatters 1:4.55). Only the critical point and the averaging radius carry the size effect; plasticity at short lives is ignored (the calibration keeps N_f > 3e4).
* **Numerics.** The authors' degree-6 polynomial fits of the averaged quantities vs r_c are used as given (R² > 99.9 %). Averaging is a plain node mean in a sphere or disc on their ANSYS meshes (0.01-0.1 mm). Step 2 (`ecp_fe.py`) recomputes these factors with our NGSolve models and `cplane`.
* **Visualisation.** Calibration and prediction are reported separately throughout. Figures: `ecp_scatter_vs_rc`, `ecp_fs_vs_n`, `ecp_life_hotspot_vs_ecp`, `ecp_series_bias`.

Data: github.com/achiocca1/ECP @ 9d51152 (GPL-3, fetched by `ecp_data.fetch()`, not redistributed).
