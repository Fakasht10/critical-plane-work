"""Valean, Foti, Berto, Marsavina, Theor. Appl. Fract. Mech. 131 (2024) 104417 (CC BY 4.0):
"Static and fatigue behavior of 3D printed smooth and notched PLA and short carbon fibers reinforced PLA".

FDM (Prusa i3 MK3), nozzle 0.4 mm, 215 C, layer 0.2 mm, raster -45/+45, 4 contours, 100 % infill, flat build.
All specimens with the same printing parameters. NOT the same batch: the smooth fatigue results are taken from
the authors' earlier study [39] (paper p. 7); static and notched fatigue are from this work.
Geometry (Fig. 1): width 20 mm, thickness 4 mm (from Delta F / Delta sigma_nom = 40 mm^2 = 10 x 4), length 91.2 mm.
  smooth fatigue: hourglass, R = 45 mm, minimum width 10 mm;
  notched: double V-notch, net width 10 mm (depth 5 mm per side), opening 2a = 30/60/90/120 deg (static),
  2a = 90 deg (fatigue); root radius not reported ("sharp").
Fatigue: R = 0.1, f = 5 Hz; Delta sigma_nom = Delta F / net area. Run-outs at 2e6 marked True.
Tensile (Table 2): E [MPa], sigma_UTS [MPa], strain at break [%].
Values transcribed from the paper's Tables 2-4.
"""
TENSILE = {"PLA": dict(E=3250.60, UTS=41.59, eb=2.28), "PLA-CF": dict(E=3203.30, UTS=33.48, eb=4.32)}

STATIC = {  # Table 3: critical load F_max [N] per opening angle 2a [deg], 5 specimens each
    "PLA": {30: [1513, 1498, 1483, 1441, 1560], 60: [1203, 1217, 1181, 1228, 1195],
            90: [1288, 1329, 1279, 1270, 1386], 120: [1711, 1691, 1705, 1742, 1704]},
    "PLA-CF": {30: [938, 921, 915, 906, 925], 60: [743, 736, 743, 734, 752],
               90: [796, 794, 852, 835, 818], 120: [1049, 1015, 1062, 1036, 1063]},
}

# Table 4: (Delta F [N], Delta sigma_nom [MPa], cycles, run-out)
FATIGUE = {
    "PLA": {
        "smooth": [(990, 24.8, 63200, False), (990, 24.8, 18543, False), (810, 20.3, 86415, False),
                   (675, 16.9, 237635, False), (576, 14.4, 648651, False), (990, 24.8, 41064, False),
                   (810, 20.3, 83950, False), (675, 16.9, 199541, False), (504, 12.6, 1987079, False),
                   (576, 14.4, 601213, False), (504, 12.6, 1046264, False), (810, 20.3, 107524, False),
                   (675, 16.9, 279654, False), (990, 24.8, 29945, False), (504, 12.6, 2000000, True)],
        "notched": [(576, 14.4, 26083, False), (576, 14.4, 30550, False), (504, 12.6, 57397, False),
                    (504, 12.6, 64355, False), (378, 9.5, 228876, False), (378, 9.5, 174834, False),
                    (306, 7.7, 514019, False), (306, 7.7, 727176, False), (576, 14.4, 22728, False),
                    (504, 12.6, 50591, False), (378, 9.5, 257756, False), (306, 7.7, 873093, False),
                    (252, 6.3, 1765445, False), (252, 6.3, 1823068, False), (252, 6.3, 2000000, True)],
    },
    "PLA-CF": {
        "smooth": [(675, 16.9, 17252, False), (675, 16.9, 12554, False), (675, 16.9, 13744, False),
                   (576, 14.4, 48415, False), (576, 14.4, 41997, False), (576, 14.4, 45868, False),
                   (405, 10.1, 393993, False), (405, 10.1, 400530, False), (405, 10.1, 572644, False),
                   (315, 7.9, 1754429, False), (315, 7.9, 2000000, True), (315, 7.9, 2000000, True),
                   (504, 12.6, 126705, False), (504, 12.6, 109671, False), (504, 12.6, 138463, False)],
        "notched": [(576, 14.4, 857, False), (504, 12.6, 1395, False), (378, 9.5, 14467, False),
                    (306, 7.7, 27509, False), (252, 6.3, 59316, False), (252, 6.3, 38634, False),
                    (252, 6.3, 37476, False), (306, 7.7, 31368, False), (306, 7.7, 24950, False),
                    (378, 9.5, 4823, False), (576, 14.4, 695, False), (576, 14.4, 578, False),
                    (504, 12.6, 872, False), (504, 12.6, 971, False), (378, 9.5, 4649, False)],
    },
}

# Paper's own results (for comparison): nominal S-N (P.S. 50 %) N = 10^a dS^-k; SED summary scatter T_sigma (5-95 %)
PAPER = {
    # Fig. 5a/b legends (re-read at 300 dpi by audit): log10 N = a - k log10 dS; T_nom = paper's scatter index per curve
    "PLA": dict(smooth=(12.11, 5.45), notched=(10.31, 5.07), R0_static=0.81, R0_fat=(0.91, 0.97), T_SED=1.24,
                T_nom=dict(smooth=1.36, notched=1.19)),
    "PLA-CF": dict(smooth=(12.14, 6.47), notched=(9.09, 5.42), R0_static=0.30, R0_fat=(0.52, 0.69), T_SED=1.29,
                   T_nom=dict(smooth=1.16, notched=1.45)),
}
W, T, NET, DEPTH, LEN, R_HOURGLASS = 20.0, 4.0, 10.0, 5.0, 91.2, 45.0
