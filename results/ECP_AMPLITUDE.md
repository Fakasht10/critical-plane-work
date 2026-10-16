# Shear-amplitude measure for 90 deg out-of-phase loading (ECP, En3B V-notched bars)

**Summary.**
* **Best measure: path length / 4.** On the 41 out-of-phase specimens (6 series), the mean prediction error E_N is:
  * MCC: -0.51 (lives 3.2x too long, unsafe).
  * MRH: -0.33.
  * path length / 4: -0.17.
  * their MPHP corrector: +0.41.
  path length / 4 is a prediction: nothing was fitted to these data. Calibration and all proportional series are unchanged (identical calibration for every measure).
* **Within a factor of 3:** 73 % with path length / 4, against 41 % with MCC and 56 % with their corrector.
* **Remaining split by R** (path length / 4): R = -1 mean -0.25, R = 0 mean -0.10.
* **Physics.** Out-of-phase loading of a ductile steel causes extra cyclic (non-proportional) hardening and more active slip systems. MCC sees only the size of the shear path, not how round it is. MRH and the path length grow with roundness (by up to sqrt(2) and pi/2 for a circular path), which is the direction the data require. A material-specific non-proportionality constant (Kanazawa / Itoh type) would be the next refinement, but it needs an extra calibration series.

## All 19 prevision series and by family (r_c = 0.2 mm)
Cells: mean E_N / std / % within factor 2 / % within factor 3.

| variant | all 19 | in-phase | 90 deg | 90 deg, R = -1 | 90 deg, R = 0 |
|---|---|---|---|---|---|
| MCC | -0.19 / 0.43 / 47 / 70 | -0.15 / 0.39 / 42 / 72 | -0.51 / 0.44 / 29 / 41 | -0.73 / 0.37 / 11 / 21 | -0.32 / 0.42 / 45 / 59 |
| MRH | -0.14 / 0.39 / 52 / 76 | -0.15 / 0.39 / 42 / 72 | -0.33 / 0.40 / 49 / 61 | -0.44 / 0.43 / 42 / 53 | -0.23 / 0.35 / 55 / 68 |
| path length / 4 | -0.09 / 0.37 / 52 / 79 | -0.15 / 0.39 / 42 / 72 | -0.17 / 0.39 / 49 / 73 | -0.25 / 0.43 / 53 / 58 | -0.10 / 0.34 / 45 / 86 |
| MCC x their MPHP corrector | +0.07 / 0.42 / 49 / 74 | -0.15 / 0.39 / 42 / 72 | +0.41 / 0.38 / 37 / 56 | +0.29 / 0.40 / 42 / 58 | +0.51 / 0.35 / 32 / 55 |
| their factors + MPHP (paper) | -0.00 / 0.44 / 45 / 73 | -0.25 / 0.40 / 35 / 62 | +0.25 / 0.49 / 29 / 61 | +0.03 / 0.55 / 21 / 58 | +0.44 / 0.35 / 36 / 64 |

## Per 90 deg series
| series | MCC | MRH | path length / 4 | MCC x their MPHP corrector | their factors + MPHP (paper) |
|---|---|---|---|---|---|
| V-notch r=0.2, 90 deg, R=-1 | -0.98 ± 0.30 | -0.81 ± 0.30 | -0.61 ± 0.30 | -0.02 ± 0.29 | -0.54 ± 0.30 |
| V-notch r=0.2, 90 deg, R=0 | -0.27 ± 0.34 | -0.27 ± 0.34 | -0.22 ± 0.34 | +0.51 ± 0.35 | +0.44 ± 0.35 |
| V-notch r=1.25, 90 deg, R=-1 | -0.45 ± 0.28 | -0.07 ± 0.28 | +0.12 ± 0.28 | +0.61 ± 0.28 | +0.55 ± 0.27 |
| V-notch r=1.25, 90 deg, R=0 | +0.02 ± 0.22 | +0.02 ± 0.22 | +0.15 ± 0.22 | +0.73 ± 0.22 | +0.67 ± 0.22 |
| V-notch r=4, 90 deg, R=-1 | -0.81 ± 0.34 | -0.51 ± 0.35 | -0.31 ± 0.35 | +0.23 ± 0.36 | -0.00 ± 0.38 |
| V-notch r=4, 90 deg, R=0 | -0.70 ± 0.32 | -0.45 ± 0.32 | -0.26 ± 0.32 | +0.28 ± 0.33 | +0.20 ± 0.33 |

## Amplitude ratio FS(measure) / FS(MCC) at r_c = 0.2 mm, mean over specimens
| series | MRH | path length / 4 | their MPHP corrector |
|---|---|---|---|
| V-notch r=0.2, 90 deg, R=-1 | 1.097 | 1.217 | 1.675 |
| V-notch r=0.2, 90 deg, R=0 | 1.000 | 1.028 | 1.517 |
| V-notch r=1.25, 90 deg, R=-1 | 1.232 | 1.363 | 1.771 |
| V-notch r=1.25, 90 deg, R=0 | 1.002 | 1.074 | 1.468 |
| V-notch r=4, 90 deg, R=-1 | 1.177 | 1.306 | 1.749 |
| V-notch r=4, 90 deg, R=0 | 1.144 | 1.270 | 1.687 |

Their corrector (1.46-1.77) exceeds pi/2 = 1.571, the largest possible length/MCC ratio on a single plane, in 4 of 6 series. So their MPHP is not simply the half-perimeter on the MCC plane; the paper does not define it fully.

## Checklist
* **Theory.** FS = (τ_a / G)(1 + k σ_n,max / σ_y) on the plane that maximises it, separately for each measure. All three measures equal the half-range on a straight (proportional) path, so only the 90 deg series can change.
* **Physics.** Linear elastic, with no plasticity correction. Non-proportional hardening enters only through the shape of the path.
* **Numerics.**
  * Path: 72 samples per cycle.
  * MRH: 48-angle scan plus golden-section polish.
  * Path length: the closed polygon.
  * Kernel: verified against closed forms (ellipse, line) and a 30k-normal brute force (`tests/test_cplane.py`).
  * Screening bounds: √2·Tresca (MRH) and the Tresca path length (length).
* **Visualisation.** `ecp_amplitude.png`:
  * left: the critical-plane shear path for one specimen, with the three amplitude radii;
  * right: the per-series prediction error for each variant.

