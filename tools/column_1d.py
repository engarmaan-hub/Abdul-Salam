"""
Independent 1-D check of the thermal and settlement model (no ANSYS).

Solves transient heat conduction with phase change (enthalpy method, the
same ENTH / KXX tables, surface functions, geothermal flux and initial
profile as ansys/01_thermal_transient.mac) in two vertical columns:

  natural  natural ground, surface = dibiao (-0.48 + 11 sin(...) + 5.94e-6 t)
  embank   7 m embankment far from the abutment, top = lumian (1.5, 12.5)

and then estimates the far-field settlement the way the ANSYS model defines
it (paper mode, Eq. 4 + load compression + creep, 1-D / oedometric).

Use it to know what ANSYS should give far from the abutment (d ~ 20 m in
Fig. 14), and to check the thermal inputs against Fig. 9 (natural field).

Usage:  python3 tools/column_1d.py [--years 15] [--spinup 60] [--plot]

--spinup N starts from an N-year spin-up of the natural ground (no warming)
instead of the ICT profile of the macro.
"""
import math
import sys

import numpy as np

# ---------------------------------------------------------------- inputs
TPTS = np.array([-20, -10.01, -10, -5.01, -5, -3.01, -3, -2.01, -2, -1.01, -1, -0.51,
                 -0.5, -0.21, -0.2, 0, 0.01, 40])
C = {1: [977.2, 977.2, 1152, 1152, 1776, 1776, 3362, 3362, 5572, 5572, 18693, 18693, 30344, 30344, 66718, 66718, 1266, 1266],
     2: [810, 810, 992.5, 992.5, 1231, 1231, 1278, 1278, 3126, 3126, 5278, 5278, 14212, 14212, 99886, 99886, 1044, 1044],
     3: [1222, 1222, 1737, 1737, 2702, 2702, 6640, 6640, 6678, 6678, 11864, 11864, 35903, 35903, 130278, 130278, 1608, 1608],
     4: [981.8, 981.8, 1476, 1476, 2364, 2364, 3658, 3658, 6160, 6160, 16080, 16080, 39562, 39562, 42466, 42466, 1272, 1272],
     5: [706.6, 706.6, 820.2, 820.2, 973.5, 973.5, 1004, 1004, 2156, 2156, 3497, 3497, 9060, 9060, 62405, 62405, 861.7, 861.7]}
RHO = {1: 1800, 2: 1900, 3: 1600, 4: 1800, 5: 2060}
KF = {1: 6552, 2: 9405, 3: 7632, 4: 6552, 5: 5040}      # J/(m h degC), frozen
KU = {1: 5760, 2: 6897, 3: 5112, 4: 5760, 5: 4140}      # unfrozen
QGEO = 72.0                                              # J/(m2 h)
ICT_Y = [-30, -16, -15, -14, -13, -12, -11, -10, -9, -8, -7, -6, -5, -4, -3, -2, -1.75, -1, -0.75,
         -0.5, -0.25, 0, 0.25, 4, 4.25, 7, 10]
ICT_T = [-0.25, -0.25, -0.27, -0.29, -0.31, -0.33, -0.36, -0.41, -0.45, -0.50, -0.55, -0.60, -0.65,
         -0.69, -0.80, -0.80, 0.20, 0.20, 0.25, 0.25, 1.0, 1.0, 2.0, 2.0, 3.0, 3.0, 3.0]
YB, DZ = -30.0, 0.25

# Table 3 (mechanical): a1 b1 a2 b2 alpha% A B C D ; mat 5 = cobbly fill
MECH = {1: (34, 30, 0.42, -0.007, 1, 2e-14, 2.06, 0.32, 0.10),
        2: (36, 30, 0.42, -0.007, 5, 2e-14, 2.06, 0.32, 0.10),
        3: (28, 26, 0.40, -0.008, 12, 2.98e-7, 1.07, 0.07, 0.52),
        4: (140, 108, 0.25, -0.004, 1, 2.86e-7, 1.07, 0.06, 0.52),
        5: (61, 53, 0.35, -0.007, 1, 3.68e-7, 1.18, 0.08, 0.54)}
GAMMA = {m: RHO[m] * 9.81 for m in RHO}                  # N/m3 (dry density, as in ANSYS)
Q_TRAF = 60.1e3                                          # Pa


def enth(m):
    """ANSYS ENTH table: trapezoidal integral of rho*C over the 18 points."""
    c, h = np.array(C[m]), [0.0]
    for i in range(1, len(TPTS)):
        h.append(h[-1] + RHO[m] * (c[i] + c[i - 1]) / 2 * (TPTS[i] - TPTS[i - 1]))
    return np.array(h)


WARM = 5.94e-6     # degC/h (0.052 degC/yr)


def surface(t, mean, amp, warming=True):
    return mean + amp * math.sin(2 * math.pi * t / 8760 + math.pi) + (WARM * t if warming else 0.0)


def material(y):
    if y > 0:
        return 5
    d = -y
    return 1 if d < 0.5 else 2 if d < 2 else 3 if d < 8 else 4


def run(top, mean, amp, years, dt=2.0, init=None, warming=True):
    """Explicit enthalpy FD. Returns y (nodes), times, T history (monthly)."""
    y = np.arange(top, YB - 1e-9, -DZ)
    n = len(y)
    # half-cell materials above/below each node
    up = [material(yy + DZ / 2) for yy in y]
    dn = [material(yy - DZ / 2) for yy in y]
    up[0] = dn[0]
    dn[-1] = up[-1]
    # combined H(T) per node (piecewise linear, as ANSYS interpolates ENTH)
    htab = {m: enth(m) for m in C}
    groups = {}
    for i in range(n):
        groups.setdefault((up[i], dn[i]), []).append(i)
    ghs = {g: 0.5 * (htab[g[0]] + htab[g[1]]) for g in groups}
    gidx = {g: np.array(v) for g, v in groups.items()}

    def T_of_H(H):
        T = np.empty(n)
        for g, idx in gidx.items():
            hs = ghs[g]
            h = H[idx]
            T[idx] = np.interp(h, hs, TPTS)
            lo, hi = h < hs[0], h > hs[-1]          # linear extrapolation
            T[idx[lo]] = TPTS[0] + (h[lo] - hs[0]) / ((hs[1] - hs[0]) / (TPTS[1] - TPTS[0]))
            T[idx[hi]] = TPTS[-1] + (h[hi] - hs[-1]) / ((hs[-1] - hs[-2]) / (TPTS[-1] - TPTS[-2]))
        return T

    def H_of_T(T):
        H = np.empty(n)
        for g, idx in gidx.items():
            H[idx] = np.interp(T[idx], TPTS, ghs[g])
        return H

    kf = np.array([KF[material(yy - DZ / 2)] for yy in y[:-1]])   # face k (cell i..i+1)
    ku = np.array([KU[material(yy - DZ / 2)] for yy in y[:-1]])
    T = np.interp(y, ICT_Y, ICT_T) if init is None else np.interp(y, init[0][::-1], init[1][::-1])
    H = H_of_T(T)
    steps_per_month = int(730 / dt)
    hist, times = [T.copy()], [0.0]
    t = 0.0
    for mth in range(years * 12):
        for _ in range(steps_per_month):
            t += dt
            T = T_of_H(H)
            T[0] = surface(t, mean, amp, warming)
            tf = 0.5 * (T[:-1] + T[1:])
            k = np.where(tf <= 0, kf, np.where(tf >= 0.01, ku, kf + (ku - kf) * tf / 0.01))
            flux = k * (T[:-1] - T[1:]) / DZ              # downward heat flux per face
            dH = np.zeros(n)
            dH[1:-1] = (flux[:-1] - flux[1:]) / DZ
            dH[-1] = (flux[-1] + QGEO) / (DZ / 2)
            H[1:] += dt * dH[1:]
        T = T_of_H(H)
        T[0] = surface(t, mean, amp, warming)
        H[0] = H_of_T(T)[0]
        hist.append(T.copy())
        times.append(t)
    return y, np.array(times), np.array(hist)


def spinup(years=60):
    """Natural ground run without warming from ICT: an initial field in
    periodic equilibrium with the surface function (autumn, t = 0 phase)."""
    y, _, T = run(0.0, -0.48, 11.0, years, warming=False)
    return y, T[-1]


def thaw_depth(y, T, below=0.0):
    """Deepest node with T >= 0 among nodes at y <= below (0 if none)."""
    idx = np.where((y <= below + 1e-9) & (T >= 0))[0]
    return -y[idx].min() if len(idx) else 0.0


def settlement_1d(y, T, thaw_ref=0, fmax=None, fill_thaw=1.0):
    """Far-field settlement components (cm) relative to the geostatic state.

    thaw_ref = 0 (paper mode): strain = -alpha wherever T > 0 (from frozen)
    thaw_ref = 1: soil thawed at t = 0 (fill, top 1.9 m of ground) starts
                  strain-free, so refreezing it gives heave (+alpha)
    compression: (fill weight + traffic) / M(T), oedometric
    returns dict: emb_thaw, emb_comp, fnd_thaw, fnd_comp
    """
    out = dict(emb_thaw=0.0, emb_comp=0.0, fnd_thaw=0.0, fnd_comp=0.0)
    top = y[0]
    for i in range(len(y) - 1):
        ym = 0.5 * (y[i] + y[i + 1])
        tm = 0.5 * (T[i] + T[i + 1])
        m = material(ym)
        a1, b1, a2, b2, alpha, *_ = MECH[m]
        ta = -tm if tm < 0 else 0.0
        E = (a1 + b1 * ta) * 1e6
        nu = a2 + b2 * ta
        M = E * (1 - nu) / ((1 + nu) * (1 - 2 * nu))      # constrained modulus
        f = min(1.0, max(0.0, tm / 0.1))                   # thawed fraction
        if fmax is not None:                               # irreversible thaw
            f = fmax[i]
        if ym > 0:
            f *= fill_thaw
        thawed_at_start = ym > -1.9
        if thaw_ref == 1 and thawed_at_start:
            f -= 1.0
        thaw = alpha / 100 * f * DZ
        if ym > 0:
            out["emb_thaw"] += 100 * thaw
            out["emb_comp"] += 100 * (GAMMA[5] * (top - ym) + Q_TRAF) / M * DZ
        else:
            out["fnd_thaw"] += 100 * thaw
            out["fnd_comp"] += 100 * (GAMMA[5] * top + Q_TRAF) / M * DZ
    return out


def main():
    years = int(sys.argv[sys.argv.index("--years") + 1]) if "--years" in sys.argv else 15
    print("1-D enthalpy model, same inputs as 01_thermal_transient.mac\n")

    init = None
    if "--spinup" in sys.argv:
        ns = int(sys.argv[sys.argv.index("--spinup") + 1])
        ys, Ts = spinup(ns)
        # ground from the spin-up, fill keeps the ICT value (placed warm)
        yy = np.concatenate([[10.0, 4.25], ys])
        tt = np.concatenate([[3.0, 3.0], Ts])
        init = (yy, tt)
        print(f"initial field: {ns}-year spin-up of the natural ground (no warming)")
        print("   y (m)   ICT    spin-up")
        for d in (1, 2, 2.5, 3, 4, 5, 7, 10, 12, 15, 20, 25, 30):
            print(f"  {-d:6.1f} {np.interp(-d, ICT_Y, ICT_T):6.2f} {np.interp(-d, ys[::-1], Ts[::-1]):7.2f}")
        print()
    yN, tN, TN = run(0.0, -0.48, 11.0, years, init=init)
    print("NATURAL GROUND (dibiao)")
    print(" year  max thaw depth (m)  permafrost-table T at 3 m in Oct")
    for yr in range(1, years + 1):
        sl = TN[(yr - 1) * 12 + 1: yr * 12 + 1]
        mx = max(thaw_depth(yN, T) for T in sl)
        i3 = np.argmin(abs(yN + 3))
        if yr in (1, 2, 3, 5, 10, 15) or yr == years:
            print(f" {yr:4d}  {mx:10.2f}  {TN[yr * 12][i3]:22.2f}")
    oct15 = TN[min(years, 15) * 12]
    print("\n October of year 15, natural ground (compare paper Fig. 9, K2)")
    fig9 = {3: (-0.8, -0.75), 5: (-0.9, -0.75), 10: (-0.65, -0.4), 15: (-0.35, -0.2), 20: (-0.3, None)}
    print("   depth   1-D model   paper numerical   paper measured   (Fig. 9, read off)")
    for d in (0, 1, 2, 3, 4, 5, 10, 15, 20, 25):
        i = np.argmin(abs(yN + d))
        pn, pm = fig9.get(d, (None, None))
        f = lambda v: f"{v:6.2f}" if v is not None else "     -"
        print(f"   {d:5.1f} {oct15[i]:10.2f} {f(pn):>16} {f(pm):>16}")

    yE, tE, TE = run(7.0, 1.5, 12.5, years, init=init)
    print("\nEMBANKMENT, FAR FROM THE ABUTMENT (lumian, 7 m fill)")
    print(" year | Oct top of frozen soil | max thaw below | settlement, October (cm)")
    print("      | (m above nat. ground) | ground (m)     | mode   total  found. = thaw + load | emb. = thaw + load")
    for yr in range(1, years + 1):
        T = TE[yr * 12]
        sl = TE[(yr - 1) * 12 + 1: yr * 12 + 1]
        froz = np.where(T < 0)[0]
        line0 = yE[froz[0]] if len(froz) else float("nan")
        dg = max(thaw_depth(yE, Tm) for Tm in sl)
        if yr in (1, 3, 5, 10, 15) or yr == years:
            for mode in (0, 1):
                r = settlement_1d(yE, T, mode)
                fnd, emb = r["fnd_thaw"] + r["fnd_comp"], r["emb_thaw"] + r["emb_comp"]
                head = f" {yr:4d} | {line0:21.2f} | {dg:14.2f} |" if mode == 0 else " " * 5 + "|" + " " * 23 + "|" + " " * 16 + "|"
                print(f"{head} {'paper' if mode == 0 else 'ref=1':>5} {fnd+emb:6.1f}  {fnd:5.1f} = {r['fnd_thaw']:4.1f} + {r['fnd_comp']:4.1f}"
                      f"     | {emb:4.1f} = {r['emb_thaw']:4.1f} + {r['emb_comp']:4.1f}")
    # variants of the thaw-strain rule, year 3 and year 15, October
    print("\n Thaw-strain rule variants (far field, October):   total / foundation / embankment (cm)")
    f_cells = lambda T: np.clip(0.5 * (T[:-1] + T[1:]) / 0.1, 0, 1)
    fmax = f_cells(TE[0])
    fm = {}
    for k in range(1, len(TE)):
        fmax = np.maximum(fmax, f_cells(TE[k]))
        fm[k] = fmax.copy()
    variants = (("A  paper: frozen ref, reversible        (THAW_REF=0, IRREV=0)", {}),
                ("B  frozen ref, irreversible             (THAW_REF=0, IRREV=1)", {"irr": 1}),
                ("C  as B, fill has no thaw strain        (IRREV=1, FILL_THAW=0)", {"irr": 1, "fill": 0.0}),
                ("D  thawed-at-start strain-free          (THAW_REF=1, IRREV=0)", {"ref": 1}))
    for name, v in variants:
        row = []
        for yr in (3, 15):
            r = settlement_1d(yE, TE[yr * 12], v.get("ref", 0),
                              fm[yr * 12] if v.get("irr") else None, v.get("fill", 1.0))
            fnd, emb = r["fnd_thaw"] + r["fnd_comp"], r["emb_thaw"] + r["emb_comp"]
            row.append(f"{fnd+emb:5.1f} / {fnd:5.1f} / {emb:4.1f}")
        print(f"   {name}  yr3 {row[0]}   yr15 {row[1]}")
    print("\n paper Fig. 14 at d = 20 m: year 3  17.5 / 15.5 / 2.0 cm,"
          "  year 15  25.5 / 23.5 / 2.0 cm")
    print(" (1-D far field has no abutment heat and no 2-D stress spreading,"
          " so treat it as a check of order of magnitude)")

    if "--plot" in sys.argv:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        fig, ax = plt.subplots(1, 2, figsize=(9, 5), facecolor="#fcfcfb")
        for a, yv, TT, title in ((ax[0], yN, TN, "Natural ground"), (ax[1], yE, TE, "Embankment far field")):
            for yr, ls in ((3, "--"), (15, "-")):
                a.plot(TT[yr * 12], yv, color="#2a78d6", ls=ls, lw=2, label=f"October, year {yr}")
            a.axvline(0, color="#52514e", lw=0.8)
            a.set_xlim(-3, 10)
            a.set_title(title, loc="left")
            a.set_xlabel("Temperature (°C)")
            a.set_ylabel("y (m)")
            a.grid(color="#e4e3df")
            a.legend(frameon=False)
            for s in ("top", "right"):
                a.spines[s].set_visible(False)
        fig.tight_layout()
        fig.savefig("column_1d.png", dpi=150)
        print("\nwrote column_1d.png")


if __name__ == "__main__":
    main()
