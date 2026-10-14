"""Ezeh & Susmel (2020) double-notched PLA plate: constants and the 2D quarter face (netgen OCC).

W 25, t 3 mm, net width 6 mm (notch depth C = 9.5 mm per side), length 86 mm; U r = 3 (blunt), U r = 1
(intermediate), V 35 deg r = 0.15 (sharp). PLA E = 3476 MPa, nu = 0.33.
"""
import numpy as np

E, NU = 3476.0, 0.33
W, T, C, L = 25.0, 3.0, 9.5, 86.0                                           # C = notch depth (net width 6 mm)


def plate_face(kind, r, alpha=35.0):
    """Quarter plate x in [0, W/2], y in [0, L/2]; notch at x = W/2 centred on y = 0; root at x = W/2 - c."""
    from netgen.occ import MoveTo, Circle
    xr = W / 2 - C                                                             # notch root on the ligament
    f = MoveTo(0, 0).LineTo(xr, 0).LineTo(W / 2, 0).LineTo(W / 2, L / 2).LineTo(0, L / 2).Close().Face()
    if kind == "U":
        f = f - MoveTo(xr + r, 0).Rectangle(W / 2 - xr - r + 1.0, r).Face() - Circle((xr + r, 0), r).Face()
    else:                                                                     # V with root radius r, half angle a/2
        a = np.radians(alpha / 2)
        xc = xr + r                                                           # root circle centre on the axis
        xt, yt = xc - r * np.sin(a), r * np.cos(a)                            # tangent point on the flank
        wedge = (MoveTo(xc, 0).LineTo(W / 2 + 1.0, 0).LineTo(W / 2 + 1.0, yt + (W / 2 + 1.0 - xt) * np.tan(a))
                 .LineTo(xt, yt).Close().Face())                              # counter-clockwise (else OCC splits)
        f = f - wedge - Circle((xc, 0), r).Face()
    return f, xr
