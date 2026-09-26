"""
Hand checks for the freeze-thaw settlement model (Chen et al., Eng. Fail.
Anal. 163 (2024) 108476, Eqs. 4-9 and Table 3).

1. Prints the temperature-dependent property tables that
   ansys/02_mechanical_settlement.mac generates, so the ANSYS input can be
   checked value by value (E, nu, c, phi, Drucker-Prager alpha / sigma_y,
   creep C1).
2. 1-D thaw-settlement estimate  S = sum(alpha_i * dh_i)  for a given drop
   of the permafrost table. Use it to sanity-check the FE settlement
   under the embankment centre.

Run:  python3 tools/settlement_check.py
"""
import math

# Table 3: a1 b1 (E, MPa) | a2 b2 (nu) | a3 b3 (c, MPa) | a4 b4 (phi, deg)
#          alpha (%) | A B C D (creep)
SOILS = {
    1: ("Gravel soil 0-0.5 m",        34, 30, 0.42, -0.007, 0.00, 0.600, 31, 4.0, 1, 2e-14, 2.06, 0.32, 0.10),
    2: ("Sandy soil 0.5-2 m",         36, 30, 0.42, -0.007, 0.00, 0.600, 31, 4.0, 5, 2e-14, 2.06, 0.32, 0.10),
    3: ("Sub-clay 2-8 m",             28, 26, 0.40, -0.008, 0.15, 0.090, 22, 8.0, 12, 2.98e-7, 1.07, 0.07, 0.52),
    4: ("Weathered mudstone 8-30 m", 140, 108, 0.25, -0.004, 0.10, 0.240, 28, 11.0, 1, 2.86e-7, 1.07, 0.06, 0.52),
    5: ("Embankment fill (cobbly)",   61, 53, 0.35, -0.007, 0.03, 0.094, 23, 9.5, 1, 3.68e-7, 1.18, 0.08, 0.54),
}

# Must match the switches in 02_mechanical_settlement.mac
M_E, C_MIN, PHI_MX, SIG_U, TIM_U, T_CLMP = 1.0, 5e3, 50.0, 1000.0, 1.0, 0.1
TPT = [-20, -15, -10, -7.5, -5, -4, -3, -2, -1.5, -1, -0.75, -0.5, -0.3, -0.2, -0.1, 0, 1, 40]


def props(p, t):
    _, a1, b1, a2, b2, a3, b3, a4, b4, _, A, B, C, D = p
    ta = -t if t < 0 else 0.0                       # |T| only in frozen soil
    E = (a1 + b1 * ta ** M_E) * 1e6                 # Eq. 6  (Pa)
    nu = a2 + b2 * ta                               # Eq. 7
    c = max(C_MIN, (a3 + b3 * ta) * 1e6)            # Eq. 8  (Pa)
    phi = min(PHI_MX, a4 + b4 * ta)                 # Eq. 9  (deg)
    tn = math.tan(math.radians(phi))
    rt = math.sqrt(9 + 12 * tn * tn)                # plane-strain DP fit to MC
    alf, sy = 3 * math.sqrt(3) * tn / rt, 3 * math.sqrt(3) * c / rt
    c1 = 0.0                                        # Eq. 5 -> ANSYS time hardening
    if t < 0:
        c1 = A * C * SIG_U ** (-B) * TIM_U ** (-C) * math.exp(D / max(ta, T_CLMP))
    return E, nu, c, phi, alf, sy, c1


def creep_strain(p, sig_pa, t_h, t):
    """Eq. 5 integrated: eps_c = A (sig/SIG_U)^B (t/TIM_U)^C exp(-D/T)."""
    A, B, C, D = p[10:14]
    if t >= 0:
        return 0.0
    return A * (sig_pa / SIG_U) ** B * (t_h / TIM_U) ** C * math.exp(D / max(-t, T_CLMP))


def thaw_settlement(top, bottom, layers):
    """1-D thaw settlement (m) of soil thawed between depths top..bottom (m, +down)."""
    s = 0.0
    for z0, z1, alpha in layers:
        dh = max(0.0, min(bottom, z1) - max(top, z0))
        s += alpha / 100 * dh
    return s


if __name__ == "__main__":
    for m, p in SOILS.items():
        print(f"\nMAT {m}: {p[0]}   (thaw coeff alpha = {p[9]} %)")
        print(f"{'T':>7} {'E MPa':>8} {'nu':>6} {'c kPa':>8} {'phi':>6} {'DP alf':>7} {'sigY kPa':>9} {'creep C1':>10}")
        for t in TPT:
            E, nu, c, phi, alf, sy, c1 = props(p, t)
            print(f"{t:7.2f} {E/1e6:8.1f} {nu:6.3f} {c/1e3:8.1f} {phi:6.1f} {alf:7.3f} {sy/1e3:9.1f} {c1:10.3e}")

    # Creep magnitude check: 15 years, 100 kPa, -0.3 degC
    print("\nCreep strain after 15 y at 100 kPa, T = -0.3 degC (stress in kPa, time in h):")
    for m, p in SOILS.items():
        print(f"  MAT {m} {p[0]:<28} eps_c = {creep_strain(p, 100e3, 15 * 8760, -0.3):.2e}")

    # Ground layers below the original surface: (top, bottom, alpha %)
    ground = [(0, 0.5, 1), (0.5, 2.0, 5), (2.0, 8.0, 12), (8.0, 30.0, 1)]
    print("\n1-D thaw settlement if the permafrost table drops from 1.9 m to:")
    for z in (2.5, 3.0, 4.0, 5.0, 6.0):
        print(f"  {z:.1f} m depth -> S = {thaw_settlement(1.9, z, ground)*1000:6.0f} mm")
