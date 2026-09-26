"""
Compare the ANSYS results with Chen et al. (2024), Fig. 14 (Middle Section).

Reads the CSV files written by ansys/03_extract_results.mac:
    fig14_profiles.csv   settlement vs. distance behind the abutment
    thaw_depth.csv       deepest thawed ground under each distance (optional)

and writes
    fig14_comparison.png   2x2 figure laid out like the paper's Fig. 14
plus a table (FE vs. paper) and a 1-D thaw-settlement check on the console.

Usage
    python3 tools/compare_fig14.py [results_dir]         # default: .
    python3 tools/compare_fig14.py --demo                # test with fake data

The paper values below were read off Fig. 14 by eye (about +/-1 cm) at
every 2 m. Exact values quoted in the text: year 15, d = 0 m:
total 73 cm, foundation 61 cm, embankment 12 cm (84 % / 16 %).
"""
import csv
import os
import sys

# ---- Fig. 14a, read off the figure (cm); d = 20, 18, ..., 0 m -------------
D_PAPER = [20, 18, 16, 14, 12, 10, 8, 6, 4, 2, 0]
PAPER = {
    3: {
        "total":      [17.5, 19.5, 21.5, 24.5, 28.0, 33.5, 38.0, 42.5, 47.0, 51.0, 57.5],
        "foundation": [15.5, 17.5, 19.5, 22.0, 24.5, 29.0, 33.0, 37.0, 41.0, 44.0, 46.5],
        "embankment": [2.0, 2.0, 2.5, 3.0, 3.5, 5.0, 5.5, 6.0, 6.5, 8.0, 11.0],
    },
    15: {
        "total":      [25.5, 28.5, 31.0, 35.0, 40.0, 48.0, 54.0, 59.0, 64.0, 68.0, 73.0],
        "foundation": [23.5, 26.0, 28.5, 32.0, 36.5, 43.5, 49.0, 53.0, 57.5, 60.0, 61.0],
        "embankment": [2.0, 2.5, 2.5, 3.0, 3.5, 5.0, 5.5, 6.0, 6.5, 8.0, 12.0],
    },
}

# Table 2 / Table 3: ground layers below the natural surface (top, bottom m, alpha %)
LAYERS = [(0.0, 0.5, 1), (0.5, 2.0, 5), (2.0, 8.0, 12), (8.0, 30.0, 1)]

# reference palette (dataviz skill), light mode, slots 1-3; text in neutral ink
C_TOTAL, C_FOUND, C_EMB = "#2a78d6", "#eb6834", "#1baf7a"
INK, INK2, GRID, SURFACE = "#0b0b0b", "#52514e", "#e4e3df", "#fcfcfb"


def read_csv(path):
    with open(path, newline="") as f:
        rows = list(csv.reader(f))
    head = [h.strip() for h in rows[0]]
    return [dict(zip(head, (float(v) for v in r))) for r in rows[1:] if r]


def fe_profile(rows, year, month=10):
    sel = sorted((r for r in rows if int(r["yr"]) == year and int(r["mon"]) == month),
                 key=lambda r: -r["d"])
    return sel


def thaw_settlement(depth):
    """1-D: sum(alpha * thawed thickness) for ground thawed from 0 to depth (m)."""
    return sum(a / 100 * max(0.0, min(depth, z1) - z0) for z0, z1, a in LAYERS)


def make_demo(path):
    """Fake FE output (paper values +/- a bias) to test the script."""
    with open(path, "w") as f:
        f.write("yr,mon,d,tot,fnd,emb,fnd_pct,emb_pct,tot_ft,fnd_ft\n")
        for yr in range(1, 6):
            k = min(1.0, 0.55 + 0.03 * yr)
            for d in range(21):
                t3 = [PAPER[15][q][D_PAPER.index(2 * round(d / 2))] for q in ("total", "foundation")]
                tot, fnd = 0.9 * k * t3[0], 0.9 * k * t3[1]
                for mon in (4, 10):
                    f.write(f"{yr},{mon},{d},{tot:.2f},{fnd:.2f},{tot-fnd:.2f},"
                            f"{fnd/tot*100:.2f},{(tot-fnd)/tot*100:.2f},{tot-8:.2f},{fnd-5:.2f}\n")


def main():
    args = [a for a in sys.argv[1:] if not a.startswith("--")]
    folder = args[0] if args else "."
    f14 = os.path.join(folder, "fig14_profiles.csv")
    if "--demo" in sys.argv:
        f14 = "demo_fig14_profiles.csv"
        make_demo(f14)
        print("DEMO MODE - fake FE data, only for testing the script\n")
    rows = read_csv(f14)
    years_fe = sorted({int(r["yr"]) for r in rows})
    last = years_fe[-1]
    YEARS = (3, last) if last != 3 else (3,)
    print(f"FE years available: 1-{last}; compared: {YEARS} "
          f"(paper values exist for years 3 and 15 only)\n")

    # ---------------- table: FE vs paper ----------------
    print("Settlement, October (cm)   FE vs paper Fig. 14")
    print(f"{'yr':>3} {'d':>4} | {'total':>13} | {'foundation':>13} | {'embankment':>13}")
    for yr in YEARS:
        fe = {int(r["d"]): r for r in fe_profile(rows, yr)}
        for d in (0, 4, 10, 20):
            i = D_PAPER.index(d)
            p = PAPER.get(yr)
            r = fe[d]
            ref = (lambda k: f"{p[k][i]:5.1f}") if p else (lambda k: "   - ")
            print(f"{yr:>3} {d:>4} | {r['tot']:5.1f} / {ref('total')} | "
                  f"{r['fnd']:5.1f} / {ref('foundation')} | {r['emb']:5.1f} / {ref('embankment')}")

    # ---------------- 1-D thaw check ----------------
    td = os.path.join(folder, "thaw_depth.csv")
    if os.path.exists(td):
        th = read_csv(td)
        print("\n1-D check: thaw-strain part of foundation settlement, October")
        print(f"{'yr':>3} {'d':>4} {'thaw depth m':>13} {'sum(alpha*h) cm':>16} {'FE foundation cm':>17}")
        for yr in YEARS:
            fe = {int(r["d"]): r for r in fe_profile(rows, yr)}
            for r in th:
                if int(r["year"]) == yr and int(r["month"]) == 10 and int(r["d_m"]) in (0, 5, 10, 20):
                    d = int(r["d_m"])
                    s = thaw_settlement(r["thaw_depth_m"]) * 100
                    print(f"{yr:>3} {d:>4} {r['thaw_depth_m']:13.2f} {s:16.1f} {fe[d]['fnd']:17.1f}")
        print("  FE foundation also holds the embankment-load compression and creep,")
        print("  so it should sit a few cm above sum(alpha*h).")

    # ---------------- figure ----------------
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        print("\nmatplotlib not installed - skipping the figure (pip install matplotlib)")
        return

    plt.rcParams.update({"font.size": 10, "axes.edgecolor": INK2, "axes.labelcolor": INK,
                         "xtick.color": INK2, "ytick.color": INK2, "text.color": INK})
    fig, ax = plt.subplots(2, 2, figsize=(11, 8.2), facecolor=SURFACE)
    series = (("total", "tot", C_TOTAL, "Total"),
              ("foundation", "fnd", C_FOUND, "Foundation ground"),
              ("embankment", "emb", C_EMB, "Embankment body"))
    for j, yr in enumerate(YEARS):
        a = ax[0, j]
        fe = fe_profile(rows, yr)
        d_fe = [r["d"] for r in fe]
        for key, col, color, label in series:
            a.plot(d_fe, [r[col] for r in fe], color=color, lw=2, marker="o", ms=4,
                   label=f"{label} - ANSYS")
            if yr in PAPER:
                a.plot(D_PAPER, PAPER[yr][key], color=color, lw=1.2, ls="--", marker="o", ms=7,
                       mfc=SURFACE, mec=color, mew=1.5, label=f"{label} - paper")
        a.set_title(f"(a) Year {yr}, October", loc="left", color=INK)
        a.set_ylabel("Settlement (cm)")
        a.set_ylim(0, 80)

    for j, (key, col, color) in enumerate((("foundation", "fnd_pct", C_FOUND),
                                            ("embankment", "emb_pct", C_EMB))):
        a = ax[1, j]
        for yr, ls in zip(YEARS, ("--", "-")):
            fe = fe_profile(rows, yr)
            a.plot([r["d"] for r in fe], [r[col] for r in fe], color=color, lw=2, ls=ls,
                   label=f"Year {yr} - ANSYS")
            if yr not in PAPER:
                continue
            p = [100 * f / t for f, t in zip(PAPER[yr][key], PAPER[yr]["total"])]
            a.plot(D_PAPER, p, color=color, lw=0, marker="^" if yr == 3 else "o", ms=7,
                   mfc=SURFACE, mec=color, mew=1.5, label=f"Year {yr} - paper")
        a.set_title(f"(b) Share of {key} settlement in total", loc="left", color=INK)
        a.set_ylabel("Share of total (%)")
        a.set_ylim(70, 100) if key == "foundation" else a.set_ylim(0, 30)

    for a in ax.flat:
        a.set_facecolor(SURFACE)
        a.set_xlim(20.5, -0.5)                       # abutment on the right, as in the paper
        a.set_xlabel("Distance from abutment (m)")
        a.grid(color=GRID, lw=0.8)
        a.set_axisbelow(True)
        for s in ("top", "right"):
            a.spines[s].set_visible(False)
        a.legend(frameon=False, fontsize=8)
    fig.suptitle("Embankment settlement, Middle Section: ANSYS vs. Chen et al. (2024) Fig. 14",
                 color=INK, fontsize=12)
    fig.tight_layout()
    out = os.path.join(folder if "--demo" not in sys.argv else ".", "fig14_comparison.png")
    fig.savefig(out, dpi=150, facecolor=SURFACE)
    print(f"\nwrote {out}")


if __name__ == "__main__":
    main()
