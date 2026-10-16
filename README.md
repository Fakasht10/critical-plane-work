# Critical plane work

Fast multiaxial-fatigue critical-plane search for finite-element stress fields, an NGSolve / CalculiX toolchain, and
three validation studies on polymer and polymer-composite fatigue data. Research code; every result comes with its
checks, and the limits found by an independent audit are stated in the reports.

## Contents

| file | what |
|---|---|
| `cplane.py` | Numba kernel. Criteria: Findley, MWCM, elastic SWT, Walker, Fatemi-Socie. Proportional constant amplitude (`search_prop`) and one general non-proportional cycle (`search_path`) with a choice of shear amplitude: minimum circumscribed circle (default), maximum rectangular hull, or path length / 4. Orientation-dependent strength (build or fibre direction), hot-spot screening, `superpose` for unit load cases. |
| `tests/test_cplane.py` | 23 tests: closed forms, MCC, brute force incl. anisotropic crease, rotation invariance, bound validity, screening |
| `cplane_bench.py` | speed / accuracy vs a naive loop and a dense 1 deg grid, 1e3-1e6 points -> `results/CPLANE_SPEED.md` |
| `ccx_io.py`, `ccx_ngsolve.py` | CalculiX C3D10 deck writer / runner / .frd reader; same mesh in NGSolve P2 and ccx, Kt anchors, critical plane on both -> `results/CCX_NGSOLVE.md` |
| `cp_tcd_pla.py` + `ezeh_data.py`, `ezeh_geometry.py` | critical plane + TCD point method, notched printed PLA (Ezeh & Susmel 2020) -> `results/CP_TCD_PLA.md`, `CP_TCD_PLA_RMAX.md` |
| `cp_tcd_valean.py` + `valean_data.py` | TCD collapse / split test, PLA and PLA-CF smooth + V-notch (Valean et al. 2024) -> `results/CP_TCD_VALEAN.md` |
| `cp_duda_cfrp.py` + `duda_data.py` | non-proportional tension-torsion CFRP tubes, 3 R x 3 phases (Duda et al. 2025) -> `results/CP_DUDA_CFRP.md` |
| `ecp_data.py`, `ecp_chiocca.py` | Effective Critical Plane (Chiocca & Frendo, IJF 2024) reimplemented on their En3B data (21 notched series, axial / torsion / in- and out-of-phase): exact reproduction, predictive vs self-fit scatter, r_c and corrector sensitivity -> `results/ECP_CHIOCCA.md` |
| `ecp_fe.py` | the same specimens rebuilt in NGSolve (axisymmetric bars, 3D plates), ECP factors with `cplane`, end-to-end calibration -> prediction on our own factors -> `results/ECP_FE.md` |
| `ecp_amp.py` | shear-amplitude measure (MCC / MRH / path length) for the 90 deg out-of-phase ECP series -> `results/ECP_AMPLITUDE.md` |
| `results/TOOLCHAIN.md` | toolchain decision, CalculiX gotchas, review of Abaqus-Fatigue-Fracture-Tools |

## Quick start

```bash
pip install -r requirements.txt
python -m pytest -q tests          # ~10-40 s (first run compiles the Numba kernel)
python cp_duda_cfrp.py             # ~1 min
python cp_tcd_valean.py            # ~3 min
python cp_tcd_pla.py; PLA_R0_COLUMN=max python cp_tcd_pla.py
python ccx_ngsolve.py              # needs ccx
python cplane_bench.py             # ~20 min
```

```python
import cplane as cp
res = cp.search_prop(Sa, Sm, cp.Findley(k=0.3))     # Sa, Sm: (P, 6) Voigt [xx, yy, zz, xy, yz, zx]
res = cp.search_path(H, cp.WalkerStress(0.5))       # H: (P, T, 6), one cycle
res["value"], res["normal"], res["tau_a"], res["sn_a"], res["sn_m"]
```

## Key numbers (see the reports for the checks)

- Kernel: 4.9e4 points/s full field (1e6 points in 20 s), max error vs brute force <= 1.2e-9 (isotropic,
  proportional) and ~1.5e-4 (non-proportional paths).
- NGSolve P2 and CalculiX C3D10 on the same mesh agree to 1e-5. CalculiX's equation solver must run
  single-threaded: multi-threaded SPOOLES (ccx 2.21) intermittently returned wrong solutions.

- ECP (Chiocca & Frendo): their MATLAB outputs reproduced to 4e-15. Our NGSolve + `cplane` factors match their
  ANSYS factors within 2 % (cylindrical bars) and 3-7 % (flat plates) at r_c = 0.2 mm. With our factors, predicted
  lives for the 19 prevision series are E_N = +0.07 ± 0.42 decades (theirs -0.00 ± 0.44), 74 % within a factor
  of 3, at the same calibrated r_c of 0.22 mm.

## Open items and limits (from the audit)

- **ECP, out-of-phase loading:** the shear-amplitude measure decides the result. On the 41 specimens at 90 deg
  (no parameter fitted) the mean log10 life error is: MCC -0.51 (3.2x unsafe), MRH -0.33, path length / 4 -0.17
  (73 % within a factor of 3), their MPHP corrector +0.41. Use `amplitude="length"` for non-proportional loading
  of ductile metals. The V-notch r = 0.2 mm series at R = -1 stays 4x unsafe with every measure except their
  corrector. In two 90 deg series our global-maximum factor exceeds theirs (1.26x, 1.09x), consistent with a
  local optimum in their search. The data are steel (En3B), not printed polymer.
- **Ezeh & Susmel PLA:** it is unverified whether Table 1 sigma_A at R = 0 is the amplitude or the maximum stress.
  Under the maximum-stress reading, plain max-stress TCD is best (rms log10 0.43) and the mean-stress conclusions
  reverse. Both readings are reported.
- **Valean:** one fatigue critical distance collapses smooth + notched about as well as the paper's SED (not
  better, once the scatter-index convention is matched). With one notch geometry L acts as a fitted Kf.
- **Duda CFRP:** load-proportional parameters are capped at |r| 0.62-0.77 by the within-set scatter. Isotropic
  Findley gains little out of sample (LOSO rms log10 N 0.47 vs 0.50 nominal). Linear CLT stresses cannot drive the
  paper's fibre-plane parameter. The same R on both channels is an assumption.
- Printed-polymer tension-torsion fatigue data do not exist publicly; multiaxial validation of printed parts needs
  new tests.

## Data sources

The ECP data (github.com/achiocca1/ECP, GPL-3) are fetched by `ecp_data.fetch()` at a pinned commit into
`data/ecp` and are not redistributed. Other values are transcribed from the published tables (CC BY 4.0 where stated in the data modules):
Ezeh & Susmel, Int. J. Fatigue 136 (2020) 105583; Valean, Foti, Berto, Marsavina, Theor. Appl. Fract. Mech. 131
(2024) 104417; Duda, Smolnicki, Lesiuk, Fatemi, Arch. Civ. Mech. Eng. 25 (2025) 212. The papers themselves are not
included.

Extracted from the `elasticity_rhmp` project (fakasht10/codex), where the earlier notch / EFFM work lives.
