"""
Check the thermal tables of ansys/01_thermal_transient.mac:

  * ENTH must be the trapezoidal integral of DENS*C over the 18 MPTEMP points
    (this reproduces the macro's ENTH values to 7 digits),
  * latent heat held by each table vs. rho_d * w * 334 kJ/kg (Table 2 water
    content); a ratio < 1 means part of the water is taken as unfrozen,
  * the mat 4 peak (-0.2..0 degC) that gives mudstone the same latent heat as
    gravel (same rho_d = 1800 and w = 15 % in Table 2) -> FIX_M4.

Usage:  python3 tools/check_enthalpy.py
"""
T = [-20, -10.01, -10, -5.01, -5, -3.01, -3, -2.01, -2, -1.01, -1, -0.51, -0.5, -0.21, -0.2, 0, 0.01, 40]
C = {1: [977.2, 977.2, 1152, 1152, 1776, 1776, 3362, 3362, 5572, 5572, 18693, 18693, 30344, 30344, 66718, 66718, 1266, 1266],
     2: [810, 810, 992.5, 992.5, 1231, 1231, 1278, 1278, 3126, 3126, 5278, 5278, 14212, 14212, 99886, 99886, 1044, 1044],
     3: [1222, 1222, 1737, 1737, 2702, 2702, 6640, 6640, 6678, 6678, 11864, 11864, 35903, 35903, 130278, 130278, 1608, 1608],
     4: [981.8, 981.8, 1476, 1476, 2364, 2364, 3658, 3658, 6160, 6160, 16080, 16080, 39562, 39562, 1267, 1267, 1272, 1272],
     5: [706.6, 706.6, 820.2, 820.2, 973.5, 973.5, 1004, 1004, 2156, 2156, 3497, 3497, 9060, 9060, 62405, 62405, 861.7, 861.7]}
RHO = {1: 1800, 2: 1900, 3: 1600, 4: 1800, 5: 2060}
W = {1: .15, 2: .10, 3: .30, 4: .15, 5: .06}
CFU = {1: (977.2, 1266), 2: (810, 1044), 3: (1222, 1608), 4: (981.8, 1272), 5: (706.6, 861.7)}
NAME = {1: "gravel", 2: "sandy", 3: "sub-clay", 4: "mudstone", 5: "fill"}
L_ICE = 334000.0


def enth(m, c=None):
    c = c or C[m]
    h = [0.0]
    for i in range(1, len(T)):
        h.append(h[-1] + RHO[m] * (c[i] + c[i - 1]) / 2 * (T[i] - T[i - 1]))
    return h


def latent(m, c=None):
    cf, cu = CFU[m]
    return enth(m, c)[-1] - RHO[m] * (cf * 20 + cu * 40)


if __name__ == "__main__":
    print("mat  soil       ENTH(40 C)    latent in table   rho*w*L      ratio")
    for m in C:
        lt, lw = latent(m), RHO[m] * W[m] * L_ICE
        print(f"{m:3d}  {NAME[m]:<9} {enth(m)[-1]:.6e}   {lt:.3e}        {lw:.3e}    {lt/lw:.2f}")
    lo, hi = 1267.0, 2e5
    for _ in range(100):
        x = (lo + hi) / 2
        c = C[4][:]
        c[14] = c[15] = x
        lo, hi = (x, hi) if latent(4, c) < latent(1) else (lo, x)
    c = C[4][:]
    c[14] = c[15] = round(x)
    print(f"\nmat 4 peak for equal latent heat with gravel: {round(x)} J/(kg C)"
          f" -> latent ratio {latent(4, c) / (RHO[4] * W[4] * L_ICE):.2f}")
    print("mat 4 ENTH line 13 with the fix:",
          ",".join(f"{h:.6e}" for h in enth(4, c)[12:]))
