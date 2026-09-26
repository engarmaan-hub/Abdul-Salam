# Freeze-thaw settlement of an embankment–bridge transition (ANSYS MAPDL)

This is a sequential thermo-mechanical analysis that reproduces Chen et al., *Engineering Failure Analysis* 163 (2024) 108476. The target is **Fig. 14**: settlement at 0–20 m behind the abutment in the Middle Section, in years 3 and 15.

| Step | File | Output |
|---|---|---|
| 1 Thermal | `ansys/01_thermal_transient.mac` | `THERMAL.db`, `THERMAL.rth` |
| 2 Mechanical | `ansys/02_mechanical_settlement.mac` | `MECH.db`, `MECH.rst` |
| 3 Extract | `ansys/03_extract_results.mac` | `thaw_depth.csv`, `temp_profiles.csv`, `fig14_profiles.csv`, `settlement_history.csv`, PNG contours |
| 4 Compare | `tools/compare_fig14.py` | `fig14_comparison.png`, FE vs. paper table, 1-D check |
| – | `tools/settlement_check.py` | property tables from Eqs. 6–9, creep magnitude, 1-D thaw settlement |

Run all three macros in one working directory, in order (`/INPUT,01_thermal_transient,mac`, then 02, then 03). Then run `python3 tools/compare_fig14.py <that directory>`.

---

## Step 1: thermal analysis (what was corrected)

| # | Problem in the original script | Correction |
|---|---|---|
| 1 | Jobname was `file` | `/FILNAME,THERMAL`, so steps 2 and 3 can find the results |
| 2 | `nsel,s,loc,y,-28` selected **no nodes** because the bottom was at y = −20, so no geothermal flux was applied at all | Flux is applied at `y = YB`, with the paper's value 0.02 W/m² = **72 J/m²·h** (was 131) |
| 3 | Model depth was 20 m | **30 m** (`YB = -30`), matching Table 2 (mudstone 8–30 m) and the ICT profile |
| 4 | The embankment ended at x = −20, only ~14 m behind the abutment, but Fig. 14 needs 20 m | `XL = -50` (~44 m behind the abutment), so the 0–20 m range is clear of the adiabatic, fixed boundary |
| 5 | Layers were assigned by node bands, patched with hard-coded element numbers (6430, 4927…) that change as soon as the mesh changes | Layers assigned by **element centroid** at the Table 2 depths: fill y > 0, gravel 0–0.5, sandy 0.5–2, sub-clay 2–8, mudstone 8–30 m |
| 6 | Mat 4 (mudstone) had no latent peak between −0.2 and 0 °C | Peak = 42 466 J/kg·°C, which gives mudstone the same latent heat as gravel (Table 2: same density 1800 and water content 15 %). ENTH is regenerated to match. Switch: `FIX_M4` |
| 7 | Crushed-rock layer was on, but the paper has no mitigation measure | `USE_CR = 0` for the paper case. Turn it on afterwards to test mitigation |
| 8 | The mechanical tables had copy errors (7.83e8, a repeated row) | Removed from the thermal run; step 2 builds all properties from Table 3 |

**Already consistent with the paper (checked):**
- Table 1 boundary temperatures: ground −0.48/11, embankment top 1.5/12.5, abutment 2/16, and warming 0.052 °C/yr = 5.94e-6 °C/h.
- Table 2 conductivity, heat capacity and density. Your ENTH tables are the exact trapezoidal integral of your C tables (verified to 7 digits).

**Clock:** t = 0 is mid-October (surface temperature at its mean and falling). So:
- **October of year N** = N·8760 h. This is when thaw depth is largest, and it is the month the paper compares in.
- **April of year N** = N·8760 − 4380 h.

**Thermal checks before running step 2** (all produced by step 3, part A):
- `temp_profiles.csv`, October of year 15: natural ground within about 0.1–0.2 °C of Fig. 9; −0.3 °C at a depth of 10 m near the abutment.
- `thaw_depth.csv`: the thaw depth next to the abutment should reach about 3 m (year 3) and about 6 m (year 15), as the paper describes.

**Not modelled:** the paper's protection cone in front of the abutment (Table 1: 1.5/15). It mainly affects the ground in front of the abutment, not the embankment behind it.

---

## Step 2: mechanical analysis

**One-way coupling.** The element type changes from PLANE55 to PLANE182 (with KEYOPT(3)=2, **plane strain**). CONTA172 is switched to UX/UY. Each month, `LDREAD,TEMP` reads that month's temperature field as a body load.

### Equations → ANSYS input

| Paper | ANSYS input |
|---|---|
| Eq. 4, ε<sub>k</sub> = α when thawed, 0 when frozen | Total thermal strain `MP,THSY`: 0 below 0 °C, −α above 0.1 °C. Vertical only (1-D thaw consolidation) |
| Eq. 5, ε<sub>c</sub> = Aσ<sup>B</sup>t<sup>C</sup>e<sup>−D/T</sup> | Time-hardening creep `TB,CREEP,,,,2` with C1 = A·C·e<sup>D/\|T\|</sup>, C2 = B, C3 = C − 1. Frozen soil only; σ in kPa, t in hours |
| Eqs. 6–7, E and ν | `MPDATA,EX / PRXY` at 18 temperatures from −20 to 40 °C |
| Eqs. 8–9, c and φ | `TB,EDP` Drucker–Prager, plane-strain fit to Mohr–Coulomb |
| §3.3 boundary conditions | Sides UX = 0, bottom UY = 0, surfaces free, 60.1 kPa traffic, and the embankment may separate from the abutment |

### Load steps (so the result is defined as in Fig. 14)

- **LS1 (geostatic):** self-weight of the natural ground only; the embankment, abutment and crushed rock are weightless. **This is the settlement reference.** The ground's own self-weight compression is therefore not counted as settlement, as in the paper ("initial ground stress includes its self-weight").
- **LS2 (end of construction):** embankment and abutment weight, 60.1 kPa traffic, and the initial temperature field.
- **LS3 onward:** monthly temperatures from `THERMAL.rth`, with creep on.

### Why the defaults are `THAW_REF = 0` and `IRREV = 0` ("paper" mode)

Fig. 14 shows embankment-body settlement of **about 2 cm far from the abutment and 11–12 cm next to it, with almost no change between year 3 and year 15**. That pattern is what Eq. 4 gives when it is applied as written:
- The strain is α whenever T > 0, measured from the frozen state, so the thawed part of the 7 m fill contributes 1 % × thawed thickness.
- Far from the abutment the fill mostly refreezes, leaving about 2 cm.
- Next to the warm abutment it stays thawed: 7 cm, plus compression.

The foundation settlement near the abutment then follows the thaw depth. Sub-clay has α = 12 %, and a thaw depth of about 6 m by year 15 gives about 61 cm (paper: 61 cm).

Other settings:
- `THAW_REF = 1`: soil that is already thawed at t = 0 starts strain-free.
- `IRREV = 1`: thaw settlement never recovers and there is no frost heave. Use this for design sensitivity rather than for reproducing the paper.

---

## Step 3: extracting results

`03_extract_results.mac` writes:

- **`fig14_profiles.csv`**: for every year, April and October, and d = 0…20 m:
  - `tot`: settlement of the embankment top (y = 7)
  - `fnd`: settlement of the embankment base (y = 0)
  - `emb` = tot − fnd
  - `fnd_pct`, `emb_pct`: shares of the total (Fig. 14b)
  - `tot_ft`, `fnd_ft`: the same, measured from the end of construction (freeze-thaw part only)
- **`thaw_depth.csv`**: the deepest ground point with T ≥ 0 under each distance (d = −10 is natural ground at x = +10).
- **`temp_profiles.csv`**: T(y) at d = 0, 5, 10, 20 m and in natural ground, October of years 3 and 15.
- **`settlement_history.csv`**: settlement at d = 0, 5, 10, 20 m at every load step, showing growth over 15 years.

**Distance** d is measured from the abutment back face at road level, x = XA − d with `XA = -6.2`. For d < 0.7 m the ground column below y = −1 is the concrete foundation, so the thaw depth there refers to soil below the footing.

## Step 4: comparing with the paper

```
python3 tools/compare_fig14.py path/to/results
```

This plots your year-3 and year-15 October profiles over the paper's Fig. 14 values (read off the figure to about ±1 cm; the text gives 73 / 61 / 12 cm at d = 0 in year 15). It also prints a comparison table and a 1-D check: Σ α·(thawed thickness) from `thaw_depth.csv` against the FE foundation settlement. Use `--demo` to test the script without ANSYS results.

### How to read a mismatch

| Symptom | Where to look |
|---|---|
| Foundation settlement too low near the abutment | Thaw depth too shallow. Check the step-1 thaw depth against the paper (~3 m in year 3, ~6 m in year 15) before touching the mechanical model |
| Right thaw depth but wrong settlement | α of sub-clay (12 %) dominates. Check that Σα·h from `compare_fig14.py` explains the FE value |
| Settlement far from the abutment much lower than 15–25 cm | Load-induced compression: E of warm frozen soil (Eq. 6, `M_E`), and the LS1 reference |
| Embankment body far too large | Try `THAW_REF = 1` (no thaw strain in fill that was never frozen) |

The paper itself reports 75 cm calculated against about 100 cm measured, so aim for the **trend and the magnitudes** (the concave shape toward the abutment, foundation share > 80 %), not an exact match.

---

## Open items / what to verify in your ANSYS version

1. **The exponent m in Eq. 6 is not given in the paper.** The script uses `M_E = 1`; also run 0.5.
2. **Creep units are not stated.** The script assumes kPa and hours. With Pa the creep strain would be > 50 %, which is unrealistic. With kPa it is about 5·10⁻⁴ over 15 years, so creep is a small part of the total.
3. **Commands, not yet run in ANSYS:**
   - `TB,EDP` with `TBTEMP` combined with `TB,CREEP`. If it is rejected, set `USE_CREEP = 0`; creep is small anyway.
   - `MP,THSY` with `MP,REFT`.
   - `MP,DENS` inside `/SOLU` (the staged self-weight).
   - `NSEL,R,TEMP` in POST1.
4. **Model size.** With XL = −50 and YB = −30 the mesh is roughly 15–20 k elements, which fits the ANSYS Student limit. The thermal run is 180 load steps and the mechanical run 182.
5. **Convergence.** If a month does not converge, raise the `NSUBST` maximum or use `MSTEP = 2` (Fig. 14 extraction needs `MSTEP` to be 1, 2, 3 or 6). If the abutment drifts (it is held only through contact), set KEYOPT(12)=5 on real set 10 (the footing pair).
