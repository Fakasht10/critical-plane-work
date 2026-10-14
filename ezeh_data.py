"""Ezeh & Susmel, Int. J. Fatigue 136 (2020) 105583, notched AM PLA: Table 1 fitted S-N curves per notch / raster
angle / R, sigma_A at N_A = 2e6 (PS 50 %), inverse slope k, scatter T_sigma.

OPEN ITEM: whether sigma_A in the R = 0 column is the net stress AMPLITUDE or the MAXIMUM stress is unverified;
cp_tcd_pla.py runs both readings (PLA_R0_COLUMN=amplitude|max).
"""
import numpy as np

import ezeh_geometry as M

NA = 2e6
NET = M.W - 2 * M.C                                                          # 6 mm
TABLE1 = {  # notch: {(theta_p, R): (k, sigma_A [MPa], T_sigma)}
    "blunt": {(45, -1): (6.3, 9.2, 1.698), (45, 0): (8.3, 7.5, 1.520), (30, -1): (4.8, 6.7, 1.299), (30, 0): (7.5, 6.9, 1.435),
              (0, -1): (5.4, 8.2, 1.702), (0, 0): (7.3, 7.3, 1.335)},
    "intermediate": {(45, -1): (4.3, 4.4, 1.451), (45, 0): (5.0, 4.3, 1.442), (30, -1): (4.8, 5.1, 1.305), (30, 0): (5.4, 4.7, 1.291),
                     (0, -1): (4.1, 4.9, 1.818), (0, 0): (4.6, 4.4, 1.252)},
    "sharp": {(45, -1): (3.7, 4.1, 1.170), (45, 0): (4.5, 3.7, 1.326), (30, -1): (4.1, 4.8, 1.749), (30, 0): (4.7, 3.9, 1.265),
              (0, -1): (3.7, 4.3, 1.477), (0, 0): (5.0, 4.5, 1.280)},
}
GEOM = {"blunt": ("U", 3.0), "intermediate": ("U", 1.0), "sharp": ("V", 0.15)}


def exp_sa(notch, R, N):
    """Measured S-N (net amplitude): raster angles averaged in log space; also min/max over angles."""
    curves = [sA * (NA / N) ** (1 / k) for (th, r), (k, sA, T) in TABLE1[notch].items() if r == R]
    lc = np.log(curves)
    return np.exp(lc.mean(0)), np.exp(lc.min(0)), np.exp(lc.max(0))
