# Freeze-thaw settlement: sequential thermo-mechanical analysis (ANSYS MAPDL)

Embankment and bridge-abutment transition on warm permafrost. The method follows
Chen et al., *Engineering Failure Analysis* 163 (2024) 108476, §3.2.2–3.3, Tables 1–3.

| File | What it does |
|---|---|
| `ansys/01_thermal_transient.mac` | Step 1: 15-year transient thermal run (your script, three edits marked `!>>> SETTLEMENT`). Writes `THERMAL.db` and `THERMAL.rth`. |
| `ansys/02_mechanical_settlement.mac` | Step 2: reads the month-end temperature fields and solves for stress, deformation and settlement. Writes `settlement_history.csv`. |
| `tools/settlement_check.py` | Prints every property table that step 2 generates, plus a 1-D hand check of thaw settlement. |

## How to run

1. Run `01_thermal_transient.mac` in batch or with `/INPUT`. Check the 0 °C isotherm and compare with Fig. 9 of the paper (October, cycle 15).
2. Run `02_mechanical_settlement.mac` in the **same working directory**.
3. Open `settlement_history.csv`. Columns P1–P7 are settlements in mm (positive = downward), measured from the end of construction:

| Col | Location |
|---|---|
| P1 | embankment top next to the abutment (x = −6.85) |
| P2–P4 | embankment top at x = −10, −15, −20 |
| P5 | top of the abutment cap |
| P6 | embankment base, original ground, x = −10 |
| P7 | natural ground in front of the abutment, x = +5 |

- **P1 − P5** is the differential settlement at the bridge head (the "bump").
- **P1 − P6** is the part that comes from the embankment itself.

## How settlement is computed

**One-way coupling.** Temperature changes the soil, but deformation does not change temperature. The thermal run is solved first. Then, for every month *k*, the structural run reads `T(x, y, t_k)` with `LDREAD,TEMP` as a body load. Both runs use the same clock (hours), so month *k* is `t = 730·k` h in both. `OUTRES,ALL,5` in step 1 already writes the last substep of every month, which is what step 2 reads.

**Element change.** `ETCHG,TTS` turns PLANE55 into PLANE182. KEYOPT(3) is then set to 2 (**plane strain**); ETCHG would otherwise leave plane stress. CONTA172 is switched to KEYOPT(1)=0 (UX, UY). The embankment can open a gap from the abutment (standard contact), with friction 0.6 / 0.5 / 0.3 as in your contact pairs.

### Equations → ANSYS input

| Paper | Meaning | ANSYS implementation |
|---|---|---|
| Eq. 4, ε<sub>k</sub> = α (0 when T < 0 °C) | thaw settlement strain | Total thermal strain `MP,THSY` = 0 below `THS_1` and −α/100 above `THS_2` (0 → 0.1 °C ramp). `MP,REFT` = −1 for soil frozen at t = 0 (mats 3, 4) and +1 for soil thawed at t = 0 (mats 1, 2, 5), so each soil starts strain-free. Vertical only by default (`VERT_ONLY=1`, 1-D thaw consolidation). |
| Eq. 5, ε<sub>c</sub> = Aσ<sup>B</sup>t<sup>C</sup>e<sup>−D/T</sup> | frozen-soil creep | Differentiated in time, this is ANSYS time-hardening creep `TB,CREEP,,,,2` with C1 = A·C·e<sup>D/\|T\|</sup>, C2 = B, C3 = C − 1, C4 = 0. The temperature term goes into a temperature-dependent C1 table, which is 0 in thawed soil and uses \|T\| ≥ 0.1 °C to avoid blow-up at 0 °C. `RATE,ON` from step 2 onward. |
| Eq. 6, E = a<sub>1</sub> + b<sub>1</sub>\|T\|<sup>m</sup> | stiffness | `MPDATA,EX` on 18 temperatures from −20 to 40 °C, with \|T\| = 0 when thawed |
| Eq. 7, ν = a<sub>2</sub> + b<sub>2</sub>\|T\| | Poisson's ratio | `MPDATA,PRXY`, same temperatures |
| Eqs. 8–9, c, φ = a + b\|T\| | strength | Extended Drucker–Prager `TB,EDP,LYFUN` fitted to Mohr–Coulomb in plane strain: α = 3√3·tanφ/√(9+12tan²φ), σ<sub>Y</sub> = 3√3·c/√(9+12tan²φ). Zero dilatancy (`LFPOT` = 0). |
| §3.3 boundaries | | Sides UX = 0, bottom UY = 0, surfaces free, self weight, 60.1 kPa on the embankment top |

The load step sequence:

1. **LS1** applies self weight and traffic at the initial ICT temperature profile, with no creep. This is the end-of-construction reference state.
2. **LS2 … LS181** apply the month-end temperature fields with creep on. Settlement = UY(t) − UY(LS1).

**Where the settlement comes from.** Ice-rich sub-clay (α = 12 %) dominates, so settlement grows roughly as 12 cm for every metre the permafrost table drops into the sub-clay. `tools/settlement_check.py` gives a 1-D estimate from the 0 °C isotherm depth. Use it to check the FE result under the embankment centre: a drop from 1.9 m to 5 m gives about 365 mm.

### Reversible or irreversible thaw

- **`IRREV=0`** (default, same as the paper): the soil follows the current temperature. Soil that thaws and does not refreeze (the permafrost under the embankment) settles permanently. The active layer heaves in winter and settles back in summer.
- **`IRREV=1`**: each soil node uses the highest temperature it has reached so far. Thaw settlement then never recovers (an upper bound). Stiffness and creep are then also evaluated at that highest temperature, which is conservative.

## What you need to decide or check

1. **m in Eq. 6 is not given in the paper.** The script uses `M_E = 1`. Run 0.5 and 1.0 to see the effect. It matters only for cold soil (|T| > 1 °C); near 0 °C the result barely changes.
2. **Creep units are not stated in the paper.** The script uses σ in kPa and t in hours (`SIG_U = 1000`, `TIM_U = 1`). That gives creep strain of about 5·10⁻⁴ over 15 years in sub-clay at −0.3 °C, which is plausible. With σ in Pa, A = 3.68e-7 would give unrealistic strain (> 50 %).
3. **Geothermal flux fix.** In your original script, `nsel,s,loc,y,-28` selected no nodes (the bottom is y = −20), so no flux was applied. This is now fixed. Also note the paper uses 0.02 W/m² (72 J/m²·h), while you used 131. The value is set by `QGEO`, so re-check the thermal validation after this change.
4. **Mat 4 heat capacity** still has no latent peak between −0.2 and 0 °C (your own CHECK note). ENTH has the same gap, so fix both together if it is a typo.
5. **The old EX tables in the thermal script** have copy errors (7.83e8, and a repeated row). Step 2 deletes them and rebuilds all properties from Table 3.
6. **Commands to verify in your ANSYS version.** I wrote this without an ANSYS licence, so the scripts have not been run.
   - `TB,EDP` with `TBTEMP` combined with `TB,CREEP` (EDP creep).
   - `MP,THSY` with `MP,REFT`.
   - `*GET,…,NODE,n,BF,TEMP` (only used when `IRREV=1`; `NTEMP` is the alternative item).

   If EDP + creep is rejected, set `USE_CREEP = 0` first to get the thaw settlement, which is the dominant part.
7. **Convergence.** If a month fails to converge, use `MSTEP = 2` or a larger `NSUBST` maximum. If the abutment drifts (it is held only through contact), try KEYOPT(12)=5 on the footing pair (real set 10).
