# ECP vs point method, FS vs SWT: our pipeline on the Chiocca & Frendo / Susmel & Taylor En3B data

**Summary.**
* **Recommended: FS (path length / 4) at a point, constant d = 0.11 mm.** rms E_N on the 19 prevision series is 0.36 decades, with 80 % within a factor of 3. The lowest rms is FS (length/4), PM with a life-dependent length (0.35), a marginal gain for one more parameter. The paper's ECP-FS with MPHP gives 0.44.
* **Averaging vs point (FS, path-length amplitude, constant length).**
  * ECP: r_c = 0.220 mm, E_N -0.09 / 0.38 / 78 %.
  * Point method: d = 0.110 mm, E_N -0.05 / 0.36 / 80 %.
* **Point method ≈ ECP, with d ≈ r_c / 2** (0.11 vs 0.22 mm). This is the classic TCD relation (point at L/2, area or volume of size ~L), so the volume average adds no information the point method lacks on this data set. The PM is slightly better on 90 deg (mean -0.11 vs -0.18) and more conservative on torsion (+0.36 vs +0.16).
* **The life-dependent length x(N) barely helps.** It lowers the calibration scatter (1.17 → 1.10) but not the prediction error (rms 0.36 → 0.35), so the constant length is preferred.
* **SWT is the wrong parameter for this steel, in both forms** (all-19 rms: ECP 0.75, PM 0.62; FS ~0.39). It is blind to torsion-dominated cracking, as the paper found.
* **Reproduction of the paper's SWT result with our FE factors (r_c = 0.25 mm):**
  * calibration T 1.36 (paper 1.32);
  * prevision self-fit T 5.01 (paper 4.55);
  * uniaxial E_N +0.30 ± 0.27 (paper +0.30 ± 0.26);
  * multiaxial E_N -0.59 ± 0.68 (paper -0.49 ± 0.65).

## Calibrated lengths and scatter
| method | x const [mm] | T cal | x(N): A [mm], B | x at 1e4 / 1e6 [mm] | T cal (x(N)) |
|---|---|---|---|---|---|
| FS (MCC), ECP | 0.220 | 1.171 | 0.240, -0.150 | 0.339 / 0.170 | 1.095 |
| FS (MCC), PM | 0.110 | 1.174 | 0.120, -0.175 | 0.180 / 0.080 | 1.095 |
| FS (length/4), ECP | 0.220 | 1.171 | 0.240, -0.150 | 0.339 / 0.170 | 1.095 |
| FS (length/4), PM | 0.110 | 1.174 | 0.120, -0.175 | 0.180 / 0.080 | 1.095 |
| SWT, ECP | 0.290 | 1.280 | 0.320, -0.125 | 0.427 / 0.240 | 1.168 |
| SWT, PM | 0.150 | 1.280 | 0.170, -0.150 | 0.240 / 0.120 | 1.167 |

## Prediction, constant length (mean E_N / std / % within factor 3)
| method | all 19 | flat axial | torsion | in-phase | 90 deg | rms all |
|---|---|---|---|---|---|---|
| FS (MCC), ECP | -0.18 / 0.44 / 71 % | -0.02 / 0.31 / 91 % | +0.16 / 0.40 / 56 % | -0.15 / 0.42 / 70 % | -0.50 / 0.44 / 49 % | 0.47 |
| FS (MCC), PM | -0.13 / 0.42 / 74 % | -0.05 / 0.29 / 91 % | +0.36 / 0.40 / 56 % | -0.08 / 0.38 / 75 % | -0.42 / 0.43 / 54 % | 0.44 |
| FS (length/4), ECP | -0.09 / 0.38 / 78 % | -0.02 / 0.31 / 91 % | +0.16 / 0.40 / 56 % | -0.15 / 0.42 / 70 % | -0.18 / 0.40 / 71 % | 0.39 |
| FS (length/4), PM | -0.05 / 0.36 / 80 % | -0.05 / 0.29 / 91 % | +0.36 / 0.40 / 56 % | -0.08 / 0.38 / 75 % | -0.11 / 0.37 / 73 % | 0.36 |
| SWT, ECP | -0.22 / 0.72 / 46 % | +0.34 / 0.24 / 75 % | -1.04 / 0.36 / 0 % | -0.07 / 0.56 / 50 % | -0.96 / 0.52 / 12 % | 0.75 |
| SWT, PM | -0.13 / 0.61 / 56 % | +0.32 / 0.23 / 79 % | -0.71 / 0.36 / 44 % | +0.02 / 0.50 / 60 % | -0.76 / 0.46 / 24 % | 0.62 |

## Prediction, life-dependent length x(N) (mean E_N / std / % within factor 3)
| method | all 19 | flat axial | torsion | in-phase | 90 deg | rms all |
|---|---|---|---|---|---|---|
| FS (MCC), ECP | -0.19 / 0.45 / 67 % | -0.01 / 0.33 / 84 % | +0.12 / 0.42 / 56 % | -0.15 / 0.41 / 78 % | -0.54 / 0.46 / 37 % | 0.49 |
| FS (MCC), PM | -0.15 / 0.43 / 69 % | -0.05 / 0.30 / 86 % | +0.37 / 0.39 / 56 % | -0.07 / 0.36 / 75 % | -0.48 / 0.46 / 41 % | 0.45 |
| FS (length/4), ECP | -0.08 / 0.38 / 77 % | -0.01 / 0.33 / 84 % | +0.12 / 0.42 / 56 % | -0.15 / 0.41 / 78 % | -0.14 / 0.40 / 71 % | 0.39 |
| FS (length/4), PM | -0.04 / 0.35 / 78 % | -0.05 / 0.30 / 86 % | +0.37 / 0.39 / 56 % | -0.07 / 0.36 / 75 % | -0.07 / 0.37 / 76 % | 0.35 |
| SWT, ECP | -0.19 / 0.73 / 44 % | +0.39 / 0.26 / 65 % | -0.99 / 0.39 / 0 % | -0.03 / 0.57 / 57 % | -0.98 / 0.49 / 12 % | 0.76 |
| SWT, PM | -0.13 / 0.65 / 52 % | +0.35 / 0.25 / 70 % | -0.75 / 0.35 / 22 % | +0.06 / 0.52 / 68 % | -0.84 / 0.47 / 20 % | 0.66 |

## Paper Table 4 (SWT, r_c = 0.25 mm): per-series self-fit T, ours vs paper
Not an FE check. Within a series every specimen has the same unit field, so SWT ∝ load² and the self-fit T depends only on the test data. It confirms the data and the T formula, nothing more. The FE-sensitive checks are the calibration T and the E_N values in the summary.

| series | ours | paper |
|---|---|---|
| Axial_Hole_8_R_-1 | 1.106 | 1.106 |
| Axial_V_notch_R_-1 | 1.225 | 1.225 |
| Axial_Hole_3_5_R_-1 | 1.101 | 1.101 |
| Axial_Hole_3_5_R_0.1 | 1.266 | 1.266 |
| Axial_Hole_8_R_0.1 | 1.188 | 1.187 |
| Axial_U_notch_R_-1 | 1.112 | 1.111 |
| Axial_U_notch_R_0.1 | 1.156 | 1.156 |
| Axial_V_notch_R_0.1 | 1.598 | 1.598 |
| Multi_V_notch_02_Ph_0_R_-1 | 1.522 | 1.522 |
| Multi_V_notch_02_Ph_0_R_0 | 1.217 | 1.217 |
| Multi_V_notch_02_Ph_90_R_-1 | 1.651 | 1.650 |
| Multi_V_notch_02_Ph_90_R_0 | 1.508 | 1.508 |
| Multi_V_notch_12_Ph_0_R_-1 | 1.179 | 1.178 |
| Multi_V_notch_12_Ph_0_R_0 | 1.122 | 1.108 |
| Multi_V_notch_12_Ph_90_R_-1 | 1.356 | 1.364 |
| Multi_V_notch_12_Ph_90_R_0 | 1.355 | 1.355 |
| Multi_V_notch_4_Ph_0_R_-1 | 1.385 | 1.385 |
| Multi_V_notch_4_Ph_0_R_0 | 1.468 | 1.468 |
| Multi_V_notch_4_Ph_90_R_-1 | 1.331 | 1.330 |
| Multi_V_notch_4_Ph_90_R_0 | 1.164 | 1.162 |
| Torsion_V_notch_02_R_-1 | 1.427 | 1.426 |

## Checklist
* **Theory.**
  * FS (extended plane; MCC or path length / 4 shear amplitude) and SWT (Socie: max normal-strain amplitude plane, σ_n,max on it), both on linear-elastic tensors.
  * ECP averages the tensor over a sphere or disc of radius r_c. PM takes the tensor at distance d along the notch bisector, at the hot-spot depth for the plates.
* **Physics.**
  * En3B is ductile; FS is the physically appropriate parameter for its shear cracking.
  * Plasticity is not modelled. Calibration uses N_f > 3e4 only, but the prediction set includes short lives.
* **Numerics.**
  * Fields and meshes as in `ECP_FE.md` (converged at the calibrated radius).
  * Path: 36 samples per cycle.
  * SWT plane: 4000-normal grid plus two cap refinements (to 0.03 deg).
  * x(N): 50 × 37 grid on (A, B). Predicted life is the first crossing of P(x(N)) = c N^m.
  * The V-notch plate is finely meshed only to 0.35 mm.
* **Visualisation.**
  * `ecp_compare_calibration`: calibration scatter vs length for each method.
  * `ecp_compare_families`: prediction error by loading family, constant vs life-dependent length.

