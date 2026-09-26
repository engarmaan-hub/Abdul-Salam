# Freeze-thaw settlement of an embankment–bridge transition (ANSYS MAPDL)

This is a sequential thermo-mechanical analysis that reproduces Chen et al., *Engineering Failure Analysis* 163 (2024) 108476. The target is **Fig. 14**: settlement at 0–20 m behind the abutment in the Middle Section, in years 3 and 15.

**Full explanation of the model** (physics, every equation and parameter, the corrections, the traffic-load literature, all checks, expected results and limits): **[docs/Model_Note.md](docs/Model_Note.md)**

## Files

| Step | File | Output |
|---|---|---|
| 1 Thermal | `ansys/01_thermal_transient.mac` | `THERMAL.db`, `THERMAL.rth` |
| 2 Mechanical | `ansys/02_mechanical_settlement.mac` | `MECH.db`, `MECH.rst` |
| 3 Extract | `ansys/03_extract_results.mac` | `thaw_depth.csv`, `temp_profiles.csv`, `fig14_profiles.csv`, `settlement_history.csv`, PNG contours |
| 4 Compare | `tools/compare_fig14.py` | `fig14_comparison.png`, FE vs. paper table, 1-D check |

Run all three macros in one folder, in order (`/INPUT,01_thermal_transient,mac`, then 02, then 03). Then run `python3 tools/compare_fig14.py <folder>`.

## Checks that run without ANSYS

```
python3 tools/check_geometry.py        # rebuilds K/L/AL: closed areas, BC lines, contact pairs, mesh
python3 tools/lint_apdl.py ansys/*.mac # block balance, name clashes (incl. across RESUME), field limits
python3 tools/check_enthalpy.py        # ENTH = integral of rho*C, latent heat vs water content
python3 tools/settlement_check.py      # E, nu, c, phi, Drucker-Prager, creep tables; 1-D thaw settlement
python3 tools/traffic_load.py          # equivalent traffic / pavement loads (JTG, BZZ-100)
python3 tools/column_1d.py --plot      # independent 1-D thermal model + far-field settlement
python3 tools/compare_fig14.py --demo  # test the comparison plot
```

## Main switches

| Macro | Switch | Default | Meaning |
|---|---|---|---|
| 01 | `XL`, `YB` | −50, −30 | model extent (m) |
| 01 | `QGEO` | 72 | geothermal flux, J/m²·h (0.02 W/m²) |
| 01 | `USE_CR` | 0 | crushed-rock cooling layer (paper case: off) |
| 01 | `FIX_M4` | 1 | restore the mudstone latent-heat peak |
| 02 | `THAW_REF`, `IRREV`, `FILL_THAW` | 0, 0, 1 | thaw-strain rule (paper mode). See note §5.2 and §9.5 |
| 02 | `LOAD_CASE` | 1 | 1 = paper 60.1 kPa; 2 = pavement 16.4 + traffic 17.9 kPa; 3 = user |
| 02 | `USE_CREEP`, `M_E` | 1, 1.0 | creep on; exponent m in Eq. 6 (not given in the paper) |
| 02 | `MSTEP` | 1 | months per structural step |

## Status

The macros have **not been run in ANSYS** (no licence here). Everything else was checked in Python (see above and note §9). On the first run, watch the commands listed in note §9.7:
- `TB,EDP` together with `TB,CREEP`;
- `MP,THSY` with `REFT`;
- `MP,DENS` inside `/SOLU`.
