"""
Traffic and pavement load on the embankment top: the numbers behind the
Q_PAVE / Q_VEH choice in ansys/02_mechanical_settlement.mac.

1. Equivalent soil column (Chinese subgrade practice, JTG D30 / JTG D60
   heavy vehicle 550 kN): q = N*Q / (B*L),  B = N*b + (N-1)*m + d
2. BZZ-100 standard axle (JTG D50-2017): 100 kN, dual wheels, p = 0.70 MPa,
   d = 21.30 cm, wheel centres 31.95 cm -> vertical stress at the top of
   the subgrade under the pavement (Boussinesq, homogeneous half space;
   stiff asphalt/base layers spread the load more, so this is an upper bound)
3. Pavement dead load for typical high-grade pavement layers
4. Effect of the chosen surcharge on the 1-D load compression of the ground

Usage:  python3 tools/traffic_load.py
"""
import math

import numpy as np


def soil_column(N, Q=550e3, L=12.8, b=1.8, m=1.3, d=0.6):
    """Uniform pressure (Pa) of N side-by-side 550 kN vehicles."""
    B = N * b + (N - 1) * m + d
    return N * Q / (B * L), B


def boussinesq_dual(z, p=0.70e6, dia=0.2130, spacing=0.3195, n=160):
    """sigma_z (Pa) at depth z below the midpoint between two circular loads."""
    a = dia / 2
    r = np.linspace(0, a, n)
    th = np.linspace(0, 2 * np.pi, 2 * n)
    R, TH = np.meshgrid(r, th)
    s = 0.0
    for xc in (-spacing / 2, spacing / 2):
        x = xc + R * np.cos(TH)
        y = R * np.sin(TH)
        rho2 = x ** 2 + y ** 2
        dA = R * (r[1] - r[0]) * (th[1] - th[0])
        s += np.sum(3 * p * z ** 3 / (2 * np.pi * (rho2 + z ** 2) ** 2.5) * dA)
    return s


def main():
    print("1) Equivalent soil column, 550 kN heavy vehicle (L = 12.8 m, b = 1.8 m, m = 1.3 m, d = 0.6 m)")
    for N in (1, 2, 3):
        q, B = soil_column(N)
        print(f"   N = {N} vehicles side by side: B = {B:4.1f} m  ->  q = {q/1e3:5.1f} kPa"
              f"   (= {q/19.6e3:4.2f} m of soil at 19.6 kN/m3)")

    print("\n2) BZZ-100 dual wheel (0.70 MPa) - vertical stress at the subgrade top")
    for h in (0.3, 0.5, 0.7, 0.9, 1.2):
        s = boussinesq_dual(h)
        print(f"   pavement thickness {h:3.1f} m  ->  sigma_z = {s/1e3:6.1f} kPa (peak, static, half-space)")
    print("   field measurements on expressways: 20-60 kPa at the subgrade top (dynamic peaks)")

    print("\n3) Pavement dead load (typical high-grade asphalt pavement)")
    layers = (("asphalt concrete", 0.18, 24.0), ("cement-stabilised base", 0.36, 22.0),
              ("granular subbase", 0.20, 21.0))
    tot = 0.0
    for name, h, g in layers:
        tot += h * g
        print(f"   {name:<24} {h:4.2f} m x {g:4.1f} kN/m3 = {h*g:5.1f} kPa")
    print(f"   total pavement dead load                  = {tot:5.1f} kPa")

    print("\n4) What the surcharge does to the far-field load compression (1-D, 7 m fill, 144 kPa)")
    fill = 2060 * 9.81 * 7
    base = None
    for q in (0.0, 12e3, 20e3, 35e3, 60.1e3):
        ratio = (fill + q) / fill
        base = base or ratio
        print(f"   q = {q/1e3:5.1f} kPa -> load on ground {(fill+q)/1e3:5.0f} kPa,"
              f" compression x{ratio:4.2f} of fill-only")
    print("   (thaw settlement does not depend on q at all, so q mainly changes the"
          " few-cm load-compression part)")


if __name__ == "__main__":
    main()
