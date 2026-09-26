# Freeze–thaw settlement of an embankment–bridge transition on permafrost

**A complete note on the model: what it does, why, and how every number was chosen and checked**

Reference study: K. Chen et al., *Engineering Failure Analysis* 163 (2024) 108476 (the "paper").
Model files: `ansys/01_thermal_transient.mac`, `ansys/02_mechanical_settlement.mac`, `ansys/03_extract_results.mac`.
Check tools: everything in `tools/`.

> **How to read this note.** Sections 1–2 give the idea. Sections 3–8 go through the model in the order ANSYS runs it. Section 9 lists every check I ran and what it showed. Section 10 says what results to expect and what the model cannot do. Sections 11–13 are reference material (switches, how to run, APDL commands).
>
> **Status.** ANSYS was not available to me (it needs a licence), so the macros have **not been run in ANSYS**. Instead I rebuilt and tested each part in Python: geometry, APDL syntax, thermal tables, a 1-D thermal model with the same inputs, the property tables, and the load calculations. Section 9.7 lists the few commands you should watch the first time you run.

---

## 1. What we are computing

The goal is the **settlement of the road embankment behind a bridge abutment** on warm permafrost (Qinghai–Tibet Plateau) over 15 years, expressed as in the paper's Fig. 14:

| Quantity | Definition |
|---|---|
| **total settlement** | downward movement of the embankment top (road surface, y = 7 m) |
| **foundation-ground settlement** | downward movement of the embankment base (original ground, y = 0) |
| **embankment-body settlement** | total − foundation (compression of the 7 m fill itself) |

These are evaluated at distances d = 0 … 20 m behind the abutment, in October (the season of deepest thaw) of year 3 and year 15.

### The analysis chain (sequential, one-way coupled)

```
01_thermal_transient.mac        02_mechanical_settlement.mac        03_extract_results.mac
  transient heat conduction  →    static stress / deformation     →   CSV + contour plots
  with ice/water phase change     (thaw strain, stiffness,             (Fig. 9, 10-12, 14 form)
  15 years, monthly results       strength and creep all driven         │
  THERMAL.db / THERMAL.rth        by the temperature field)             ▼
                                  MECH.db / MECH.rst                tools/compare_fig14.py
                                                                      overlay on paper's Fig. 14
```

"**One-way**" means temperature changes the soil (it thaws, softens, creeps), but deformation does not change the temperature. This is standard for thaw-settlement studies: deformations are centimetres and do not change heat flow noticeably. It lets the two analyses run one after the other on the same mesh.

---

## 2. The physics in one page

Warm permafrost (−0.3 to −1 °C) is soil held together by ice. Settlement under a road comes from four mechanisms (paper §3.2.2 and the literature on the Qinghai–Tibet Highway):

1. **Thaw settlement.** When ice-rich soil thaws, its ice turns into water, which drains and consolidates, so the soil loses volume. The paper represents this with the **thaw settlement coefficient α**: a strain that appears when the soil thaws (Eq. 4). It is by far the largest term: sub-clay with α = 12 % loses 12 cm of height per metre thawed.
2. **Compression under load.** Embankment weight (7 m × 20.2 kN/m³ ≈ 141 kPa) plus pavement and traffic compress the ground. The stiffness depends on temperature: frozen soil is stiff, warm or thawed soil is soft (Eqs. 6–7).
3. **Creep of frozen soil.** Ice creeps under sustained stress, faster as the soil approaches 0 °C (Eq. 5).
4. **Plastic yielding.** Where stress exceeds strength (cohesion c and friction φ, which both fall as the soil warms, Eqs. 8–9), the soil yields.

**Why the abutment makes it worse.** The concrete abutment and foundation conduct heat much better than soil (10 584 vs 5 000–9 400 J/m·h·°C). The abutment surface is also warm (mean +2 °C, amplitude 16 °C). So heat flows down along the concrete into the ground behind it, the permafrost table there sinks, and a **thawed interlayer** forms that no longer refreezes in winter. Thaw settlement is therefore largest next to the abutment. That is the "bump at the bridge end" that Fig. 14 shows: 73 cm at the abutment vs 25 cm at 20 m in year 15.

---

## 3. Geometry, materials and mesh

### 3.1 Coordinates (m)

- x runs along the road. The bridge is on the right (+x) and the embankment on the left (−x).
- y is vertical. y = 0 is the original ground surface and y = 7 is the road surface.

| Part | Where |
|---|---|
| Abutment wall | x = −1.5 … 0, y = 0 … 6. The bridge seat at y = 6 (x = −1.5 … 0) is exposed |
| Abutment cap / backwall | x = −6.2 … −1.5, y = 6 … 7; back face battered from (−6.2, 6) to (−5.37, 0) |
| Foundation (footing) | x = −6.85 … 0.75, y = −7 … −1 |
| Embankment fill | y = 0 … 7, from x = XL = −50 to the abutment back face |
| Natural ground in front | y ≤ 0, x = 0 … 15 |
| Model bottom | y = YB = −30 |

**Distance behind the abutment:** d = XA − x with XA = −6.2 (the back face at road level).

### 3.2 Soil layers (Table 2 of the paper), assigned by element centroid

| Mat | Soil | Depth below original ground | Density (kg/m³) | Water content |
|---|---|---|---|---|
| 5 | embankment fill (Table 3 "cobbly soil") | above ground, 0–7 m | 2060 | 6 % |
| 1 | gravel soil | 0–0.5 m | 1800 | 15 % |
| 2 | sandy soil | 0.5–2 m | 1900 | 10 % |
| 3 | sub-clay | 2–8 m | 1600 | 30 % |
| 4 | weathered mudstone | 8–30 m | 1800 | 15 % |
| 7 | concrete (abutment, foundation) | – | 2500 | – |
| 11 | crushed rock (optional cooling layer, off for the paper case) | y = 4.5–6 | 1600 | – |
| 8, 9, 10 | contact pairs (friction 0.6 / 0.5 / 0.3) | – | – | – |

### 3.3 Two separate meshes joined by contact

The concrete (areas 1–9) and the soil (areas 10–25) have **separate nodes** at the same coordinates. They are connected by three contact pairs (CONTA172 on soil, TARGE169 on concrete):

| Pair | Concrete lines (target) | Soil lines (contact) | Friction μ | Cohesion | Mechanical behaviour (`ABUT_BOND = 1`) |
|---|---|---|---|---|---|
| behind the abutment (real set 8) | 3–7 | 43–47 | 0.6 | 85 kPa | standard: can slide and open a gap (paper) |
| in front (real set 9) | 13, 22, 23 | 52–54 | 0.5 | 85 kPa | **bonded** (element type 4) |
| under the footing (real set 10) | 8, 24–26 | 48–51 | 0.3 | 200 kPa | **bonded** (element type 4) |

- **In the thermal run**, contact passes heat through the thermal contact conductance TCC = 105 840 J/(m²·h·°C) = 29.4 W/(m²·K).
- **In the mechanical run**:
  - the behind pair lets the embankment slide on and separate from the abutment, as the paper assumes ("gaps between the embankment and the abutment are allowed");
  - the footing and front pairs are **bonded**, because the concrete has no other support (Section 12.1);
  - the contact stiffness factor is reset from the thermal 0.01 to `FKN_M = 1`.

### 3.4 Mesh

- Mapped quadrilaterals: 0.25 m in the upper 14 m and near the abutment, 0.5 m elsewhere.
- About **14 000 elements / 15 000 nodes** with XL = −50 and YB = −30 (the original 20 × 20 m model was about 6 300). This is well inside the ANSYS Student limit.
- `LESIZE` uses KYNDIV = 1 (soft divisions), so the mesher may adjust a few non-matching pairs of opposite sides. Those pairs exist in your original model too (it meshed). The extension adds none (Section 9.1).

---

## 4. Step 1: transient thermal analysis

### 4.1 Governing equation

2-D heat conduction with phase change (element PLANE55, time in **hours**):

$$\rho C(T)\,\frac{\partial T}{\partial t} = \nabla\cdot\big(k(T)\,\nabla T\big)$$

Latent heat of the pore ice is not a separate term. It is built into the **apparent heat capacity** C(T): a very large C in the phase-change range −0.5 … 0 °C holds the latent heat. For example, sub-clay has C = 130 278 J/kg·°C between −0.2 and 0 °C, against 1 608 when thawed.

**Why the ENTH (enthalpy) tables.** With only C(T), a time step that jumps across the narrow phase-change interval can skip the latent heat. ANSYS then melts ice "for free" and thaw depth is overestimated. With enthalpy $H(T)=\int \rho C\,dT$ given, ANSYS uses $\Delta H/\Delta T$ over each step, so latent heat is always conserved. Your ENTH tables are exactly the trapezoidal integral of ρ·C over the 18 table points. I verified all five materials to 7 digits (`tools/check_enthalpy.py`).

**How much latent heat each table holds** compared with all the water freezing ($\rho_d\,w\,L$, L = 334 kJ/kg):

| Mat | Latent in table (J/m³) | ρ·w·L | Ratio | Meaning |
|---|---|---|---|---|
| 1 gravel | 7.37e7 | 9.02e7 | 0.82 | ≈ 3 % unfrozen water |
| 2 sandy | 6.00e7 | 6.35e7 | 0.95 | nearly all water freezes |
| 3 sub-clay | 9.48e7 | 1.60e8 | 0.59 | clay keeps about 12 % unfrozen water: physically reasonable |
| 4 mudstone (original) | 5.81e7 | 9.02e7 | 0.64 | **missing peak at −0.2…0 °C** |
| 4 mudstone (FIX 6) | 7.37e7 | 9.02e7 | 0.82 | same as gravel, which has the same ρ and w in Table 2 |
| 5 fill | 4.06e7 | 4.13e7 | 0.98 | nearly all water freezes |

**Thermal conductivity:** frozen value up to 0 °C, unfrozen from 0.01 °C (Table 2). Frozen soil conducts better because ice conducts about 4× better than water. This "thermal offset" is one reason permafrost survives under a mean surface temperature slightly above 0 °C.

### 4.2 Boundary conditions

**Surface temperatures** (paper Eq. 10):

$$T(t) = T_0 + A\sin\!\Big(\frac{2\pi t}{8760}+\pi\Big) + 5.94\times10^{-6}\,t\qquad (t \text{ in h})$$

The last term is the paper's warming of 0.052 °C/year (2.6 °C in 50 years): 0.052/8760 = 5.94e-6 °C/h.

| Surface (Table 1) | T₀ (°C) | A (°C) | Lines | APDL table |
|---|---|---|---|---|
| natural ground | −0.48 | 11 | 27, 28 | `dibiao` |
| embankment top | 1.5 | 12.5 | 55, 56 | `lumian` |
| abutment concrete | 2 | 16 | 1, 2, 11, 12 | `qiaotai` |
| embankment slope (−0.18 / 13.5), protection cone (1.5 / 15) | – | – | not in a longitudinal section | – |

- **Bottom (y = −30):** geothermal heat flux 0.02 W/m² = **72 J/(m²·h)** (paper §3.3).
- **Sides (x = −50 and x = 15):** adiabatic (no heat flow). This is the standard far-field assumption, and why the sides must be far away.

**Clock.** At t = 0 the surface is at its annual mean and cooling, which is autumn. The coldest point is 3 months later (t = 2190 h, winter) and the warmest 9 months later (t = 6570 h, summer). So:
- the end of each model year, **t = N·8760 h, is "October"**, when the thaw depth is greatest. This is the month the paper validates (20 October) and reports.
- **"April" is t = N·8760 − 4380 h**.

### 4.3 Initial condition

`ICT` is a measured-style natural ground profile, interpolated at every node:
- +1 °C at the surface;
- +0.2 °C at 1.75 m;
- −0.8 °C at 2–3 m (just below the permafrost table at about 1.9 m);
- warming slowly with depth to −0.25 °C at 30 m;
- the new fill starts at +2 to +3 °C.

It is an autumn profile, consistent with t = 0 being autumn.

### 4.4 Solution controls (why they are there)

| Command | Purpose |
|---|---|
| `ANTYPE,TRANS` | transient heat transfer |
| `TIME, IM*730` in a `*DO` loop | one load step per month (730 h), 180 steps for 15 years |
| `DELTIM,1,0.01,24` + `AUTOTS` | 1 h start, 0.01–24 h automatic; small steps where phase change is active |
| `TINTP,0.009,,,1,0.5,-1` | time integration: Crank–Nicolson type (θ = 0.5) with oscillation control |
| `CNVTOL,TEMP,,,2` / `NEQIT,100` / `LNSRCH` / `PRED` | robust Newton iterations for the strongly nonlinear phase change |
| `OUTRES,ALL,5` | every 5th substep **and the last substep of every month**; step 2 reads the month ends |
| `USE_CR` block | crushed rock switches its k by season (convection in winter); off for the paper case |

### 4.5 What was corrected in the thermal script

| Fix | Problem | Why it matters |
|---|---|---|
| 1 | Jobname `file` | Steps 2 and 3 need a fixed name (`THERMAL`) |
| 2 | `nsel,s,loc,y,-28` selected **no nodes** (bottom was at −20), so no geothermal flux was applied; the value was also 131 instead of the paper's 72 J/m²h | Without the flux the deep ground cools slowly over 15 years, so permafrost degradation is underestimated |
| 3 | Depth 20 m → **30 m** | Table 2 layers go to 30 m; a shallow bottom distorts the deep temperature trend |
| 4 | Embankment only ~14 m behind the abutment → **44 m** (XL = −50) | Fig. 14 needs 20 m, and the adiabatic, fixed boundary must not sit inside that range |
| 5 | Layers by node bands plus hard-coded element numbers (6430, 4927…) | Element numbers change whenever the mesh changes, so wrong elements would get wrong soil. Centroid-based assignment is mesh-independent |
| 6 | Mat 4 latent peak missing | See the table in 4.1 |
| 7 | Crushed rock ON → OFF | The paper studies the unprotected transition; test mitigation separately |
| 8 | Old mechanical tables (with typos 7.83e8 and a repeated row) removed | Step 2 builds all mechanical properties from Table 3 |

---

## 5. Step 2: mechanical analysis

### 5.1 From thermal to structural model

- **`RESUME,THERMAL,db`**: same mesh, materials and contact.
- **`ETCHG,TTS`**: PLANE55 becomes **PLANE182**, then `KEYOPT(3)=2` sets **plane strain**. ETCHG would leave plane stress, which is wrong for a long embankment.
- **CONTA172** KEYOPT(1)=0 switches it to UX/UY. KEYOPT(9)=1 ignores tiny initial gaps, and KEYOPT(10)=2 updates contact stiffness each iteration.
- **Each month**, `LDREAD,TEMP,,,t,,THERMAL,rth` loads that month's nodal temperatures as a **body load**. They drive thermal strain and every temperature-dependent property.

Units: N, m, Pa, kg, and time in **hours** (the same clock as the thermal run).

### 5.2 Constitutive model, equation by equation

Table 3 constants are stored in the array `MPT(14,5)`, one column per soil. The macro builds temperature tables at 18 temperatures (−20 … +40 °C) in loops. **|T| acts only in frozen soil**: for T ≥ 0 every property equals its thawed value (a₁, a₂, a₃, a₄).

**Eq. 6, stiffness:** $E = a_1 + b_1|T|^m$ (MPa → ×10⁶ Pa). The paper does not give *m*; `M_E = 1` (linear).
**Eq. 7, Poisson's ratio:** $\nu = a_2 + b_2|T|$.

| Soil | E (MPa) at +1 / −0.3 / −1 / −5 °C | ν at +1 / −1 °C |
|---|---|---|
| gravel | 34 / 43 / 64 / 184 | 0.420 / 0.413 |
| sandy | 36 / 45 / 66 / 186 | 0.420 / 0.413 |
| sub-clay | 28 / 36 / 54 / 158 | 0.400 / 0.392 |
| mudstone | 140 / 172 / 248 / 680 | 0.250 / 0.246 |
| fill | 61 / 77 / 114 / 326 | 0.350 / 0.343 |

**Eqs. 8–9, strength:** $c = a_3 + b_3|T|$ and $\varphi = a_4 + b_4|T|$. ANSYS's Extended Drucker–Prager (`TB,EDP`) needs the pressure sensitivity α and yield stress σ_Y. They come from matching Drucker–Prager to Mohr–Coulomb in **plane strain**:

- Classic Drucker–Prager: $\sqrt{J_2} + a I_1 = k$, with $a = \dfrac{\tan\varphi}{\sqrt{9+12\tan^2\varphi}}$ and $k = \dfrac{3c}{\sqrt{9+12\tan^2\varphi}}$.
- ANSYS EDP (linear): $q + \alpha\sigma_m - \sigma_Y = 0$, with $q=\sqrt{3J_2}$ and $\sigma_m = I_1/3$.
- Multiplying the first form by √3 gives

$$\alpha = \frac{3\sqrt3\tan\varphi}{\sqrt{9+12\tan^2\varphi}},\qquad \sigma_Y = \frac{3\sqrt3\,c}{\sqrt{9+12\tan^2\varphi}}$$

Two numerical safeguards (not in the paper):
- **c ≥ `C_MIN` = 5 kPa.** Thawed gravel and sand have c = 0, which makes the free surface unstable numerically.
- **φ ≤ `PHI_MX` = 50°.** Eq. 9 gives, for example, 31 + 4·10 = 71° at −10 °C, which is unphysical. Near 0 °C, where most of the soil is, the cap has no effect.

The flow potential (`LFPOT`) is 0: no dilatancy (non-associated flow), the usual choice for soils.

**Eq. 5, creep of frozen soil:** $\varepsilon_c = A\,\sigma^B\,t^C\,e^{-D/T}$.

Differentiating in time gives the rate:

$$\dot\varepsilon_c = A\,C\,\sigma^B\,t^{C-1}\,e^{-D/T}$$

This is exactly ANSYS's implicit time-hardening law `TB,CREEP,,,,2`: $\dot\varepsilon = C_1\sigma^{C_2}t^{C_3}e^{-C_4/T}$. So:

$$C_1 = A\,C\;\text{SIG\_U}^{-B}\;\text{TIM\_U}^{-C}\;e^{D/|T|},\quad C_2 = B,\quad C_3 = C-1,\quad C_4 = 0$$

- The factor $e^{-D/T}$ is **folded into C₁(T)**, where T is °C and negative in frozen soil, so $e^{-D/T} = e^{D/|T|}$.
  - This lets C₁ = 0 in thawed soil (creep of frozen soil only).
  - |T| is held ≥ 0.1 °C (`T_CLMP`) to avoid the blow-up at 0 °C.
- **Units are not stated in the paper.** With σ in **kPa** (`SIG_U = 1000`) and t in hours, sub-clay at −0.3 °C and 100 kPa creeps about 5·10⁻⁴ in 15 years, which is plausible. With σ in Pa it would be > 50 %, which is impossible. So kPa is used and creep is a small part of the settlement.

**Eq. 4, thaw settlement strain:** $\varepsilon_k = \alpha$ when thawed, 0 when frozen.

It is entered as a temperature-dependent total thermal strain (`MP,THSY`), vertical only (`THSX = THSZ = 0`, `VERT_ONLY = 1`) because thaw consolidation is one-dimensional:

$$\text{THSY}(T) = \begin{cases}0 & T \le 0\\ -\alpha\,T/0.1 & 0<T<0.1\\ -\alpha & T\ge 0.1\end{cases},\qquad \varepsilon_{th} = \text{THSY}(T) - \text{THSY}(T_{\text{ref}})$$

The 0.1 °C ramp only helps convergence. α (%): gravel 1, sandy 5, **sub-clay 12**, mudstone 1, fill 1.

**Which state is "zero strain"** is set by the reference temperature T_ref (`MP,REFT`), through two switches:

| Switch | Meaning |
|---|---|
| `THAW_REF = 0` (paper) | every soil referenced to the frozen state (T_ref = −1 °C): wherever soil is thawed it carries −α |
| `THAW_REF = 1` | soil already thawed at t = 0 (fill, top 1.9 m of ground) starts strain-free; refreezing it gives heave (+α) |
| `IRREV = 0` (paper) | strain follows the current temperature (thaw settles, refreezing heaves back) |
| `IRREV = 1` | each soil node keeps the highest temperature reached, so thaw strain never recovers (thaw consolidation is irreversible; frost heave under a 200 kPa embankment is suppressed). Properties then also use that highest temperature, which is conservative |
| `FILL_THAW = 1` (paper) / 0 | fill gets α = 1 % / fill placed unfrozen has no ice, so no thaw strain |

Section 9.5 shows what each combination gives.

**Concrete** is linear elastic: E = 25 GPa, ν = 0.167, α_T = 1·10⁻⁵ /°C (T_ref = 2 °C). **Crushed rock** (if used) is elastic with E = 100 MPa and ν = 0.3.

### 5.3 Boundary conditions (paper §3.3)

- Sides x = XL and x = 15: **UX = 0**.
- Bottom y = YB: **UY = 0**.
- Surfaces are free.
- Gravity: `ACEL,0,9.81` (N/m³ = ρ·g).

---

## 6. Load steps and the reference state

| LS | Time (h) | Applied | Purpose |
|---|---|---|---|
| 1 | 0.5 | gravity on the **natural ground only** (fill, concrete and crushed rock made weightless with `MP,DENS` ≈ 0); soil at a frozen reference temperature | **geostatic stress**, the reference for settlement |
| 2 | 0.75 | full densities, pavement + traffic surcharge (temperatures unchanged) | construction loads |
| 3 | 1 | initial temperature profile ICT (thawed active layer and fill take thaw strain and thawed properties) | end of construction |
| 4 … 183 | 730·k | monthly temperature field from THERMAL.rth; creep on (`RATE,ON`) | 15 years of service |

Loads (LS2) and the first temperature field (LS3) are applied in **separate steps**. Each is a large change on its own, and splitting them makes both converge more easily.

**Why measure from LS1.** The paper states that "the initial ground stress includes self-weight". Real ground has already consolidated under its own weight, so that compression is not settlement. Measuring from LS1 therefore excludes it and includes:
1. the compression caused by the embankment and traffic;
2. thaw strain;
3. creep and plastic strain.

`03` also writes a "freeze-thaw only" version measured from LS3, the end of construction (`tot_ft`, `fnd_ft`).

---

## 7. Traffic and pavement load: what value to use

### 7.1 What the paper uses

"The dynamic load is equivalent to a static load with a strength of **60.1 kPa**" (paper §3.3, their ref. [39]), applied uniformly on the road surface.

### 7.2 What the literature and codes give

The calculated values come from `tools/traffic_load.py`.

| Source / method | Value | What it represents |
|---|---|---|
| **Equivalent soil column** (Chinese subgrade practice, JTG D30 with the JTG D60 heavy vehicle of 550 kN): q = N·Q/(B·L), B = N·b + (N−1)·m + d, with L = 12.8 m, b = 1.8 m, m = 1.3 m, d = 0.6 m | **17.9 kPa** (1 lane), 15.6 (2), 15.0 (3); about 0.8–0.9 m of soil | static traffic spread over the vehicle footprint, for embankment and subgrade design |
| **FHWA / AASHTO** live-load surcharge | **12 kPa** (250 psf) | highway traffic surcharge on embankments and walls |
| **Eurocode LM1-based equivalent loads** for embankments (plane strain) | **20 kPa** (embankment < 4 m), **15 kPa** (> 4 m, dual carriageway); 31 kPa heavy-vehicle strip + 9 kPa elsewhere | equivalent characteristic traffic load for geotechnical analyses |
| **BZZ-100** standard axle (JTG D50-2017): 100 kN, dual wheels, 0.70 MPa, d = 21.30 cm, 31.95 cm apart; Boussinesq stress at the subgrade top | 138 / **73** / **42** / 27 kPa under 0.3 / 0.5 / 0.7 / 0.9 m of pavement | **peak** stress under one wheel pair, very local, upper bound (stiff layers spread more) |
| Field measurements, expressway subgrade top | **20–60 kPa** (30–50 kPa typical), dropping to 10–20 kPa at 3 m depth | short **dynamic** peaks |
| **Pavement dead load** (0.18 m asphalt ×24 + 0.36 m cement-stabilised base ×22 + 0.20 m subbase ×21 kN/m³) | **16.4 kPa** | permanent weight of the pavement layers, which the 7 m fill in the model does not include |

### 7.3 Interpretation and recommendation

- **The paper's 60.1 kPa** equals the peak stress under one BZZ-100 wheel pair through about 0.55 m of pavement, and sits at the top of the measured 20–60 kPa dynamic range. Applied uniformly and permanently over the whole road, it is **conservative** (roughly 2× a realistic sustained load).
- **For long-term settlement** (thaw, consolidation, creep) the relevant load is the **sustained, spread** load: pavement dead load + equivalent traffic = 16.4 + 17.9 ≈ **34 kPa**. This is `LOAD_CASE = 2`. Repeated dynamic loading causes cumulative plastic strain, a separate mechanism this model does not include.
- **How much it matters:** the traffic load does **not** change thaw settlement, which is purely temperature-driven. It changes only load compression. Far from the abutment, going from 34 to 60 kPa raises the ground load from 176 to 202 kPa, about +15 % of a roughly 4 cm compression term, so **under 1 cm**. Next to the abutment, where the ground thaws and softens, the effect is larger but still secondary to α·(thaw depth).

| `LOAD_CASE` | Q_PAVE | Q_VEH | Total | Use |
|---|---|---|---|---|
| 1 (default) | 0 | 60.1 kPa | 60.1 kPa | reproduce the paper |
| 2 | 16.4 kPa | 17.9 kPa | 34.3 kPa | realistic sustained load (design) |
| 3 | your values | your values | – | sensitivity |

`LOAD_CAP = 1` also puts the surcharge on the road over the abutment cap (line 2), because the pavement and traffic continue over the abutment.

---

## 8. Step 3: extracting the results

`03_extract_results.mac` has two parts.

**Part A: thermal (`RESUME THERMAL.db`, `FILE THERMAL.rth`)**
- `thaw_depth.csv`: for every year, April and October, and d = 0…20 m (and natural ground at x = +10, coded d = −10), the deepest soil node with T ≥ 0 below the original ground. The paper describes about 3 m next to the abutment in year 3 and about 6 m in year 15.
- `temp_profiles.csv`: T(y) every 0.5 m at d = 0, 5, 10, 20 m and in natural ground, October of years 3 and 15 (compare Fig. 9).
- PNG contours of temperature (−2 … +2 °C bands) for October year 3, October year 15 and April year 15.

**Part B: settlement (`RESUME MECH.db`, `FILE MECH.rst`)**
- Nodes: top `NODE(XA−d, 7)` and base `NODE(XA−d, 0)` among soil nodes only (the concrete has nodes at the same coordinates).
- `fig14_profiles.csv`, columns:

| Column | Meaning |
|---|---|
| `yr, mon, d` | year, month (4 or 10), distance (m) |
| `tot, fnd, emb` | total, foundation, embankment settlement (cm), from LS1 |
| `fnd_pct, emb_pct` | shares of the total (%), as in Fig. 14b |
| `tot_ft, fnd_ft` | the same, from LS3 (freeze–thaw part only) |

- `settlement_history.csv`: every load step, total and foundation settlement at d = 0, 5, 10, 20 m.
- PNG contours of vertical displacement relative to LS1 (`LCDEF`/`LCOPER,SUB`).

Then run `python3 tools/compare_fig14.py <folder>`. It:
- plots your year-3 and year-15 profiles over the paper's Fig. 14 (values read off the figure to about ±1 cm; the text gives 73 / 61 / 12 cm at d = 0 in year 15);
- prints an FE-vs-paper table;
- prints a 1-D check, Σα·(thawed thickness), against the FE foundation settlement.

---

## 9. Checks I performed (and what they showed)

### 9.1 Geometry rebuilt in Python (`tools/check_geometry.py`)

The script replays every `K`, `L` (including `*DO` loops; ANSYS does not create a second line between the same two keypoints) and `AL` command.

- 50 keypoints, **73 lines**, 25 areas, **every area a closed 4-line loop**.
- BC lines are what the comments say:
  - 27–28 are the ground surface (x = 0 … 15);
  - 55–56 are the embankment top (x = −50 … −6.2);
  - 1, 2, 11, 12 are the exposed concrete faces.
- All three contact pairs are **geometrically coincident** (target and contact lines lie exactly on top of each other).
- **Every exterior line has a boundary condition:** surface temperature, the bottom heat flux, or an adiabatic side.
- Mapped mesh: about 14 000 elements. The extension adds no new division mismatch (Section 3.4).

### 9.2 APDL lint (`tools/lint_apdl.py`)

- Checks `*DO`/`*IF` block balance, line and field limits, and names clashing with APDL functions.
- It also catches **cross-macro clashes**, because `RESUME` restores parameters.
- **Bug found and fixed:** `03` created a component `SOILN` while an array parameter `SOILN` from `02` was restored with `MECH.db`. The component is now `SOILNODE`.
- Also made more robust: scalar variables instead of array/function arguments in `*IF`, and `DEL` renamed to `DELA`. All three macros now pass.

### 9.3 Thermal tables (`tools/check_enthalpy.py`)

- ENTH = ∫ρC dT, reproduced to 7 digits.
- Latent-heat ratios as in Section 4.1.
- The mat-4 fix value 42 466 J/kg·°C, with the regenerated ENTH line identical to the macro.

### 9.4 Independent 1-D thermal model (`tools/column_1d.py`)

- Explicit enthalpy finite differences (dz = 0.25 m, dt = 2 h).
- Uses the **same** ENTH and k tables, surface functions, 72 J/m²·h flux and ICT profile as the macro.

**Natural ground** (compare the paper's Fig. 9, natural field K2; paper values read off the figure):

| Depth (m) | 1-D, start from ICT | 1-D, start from 60-yr spin-up | Paper, numerical | Paper, measured |
|---|---|---|---|---|
| 3 | −0.43 | −0.47 | ≈ −0.80 | ≈ −0.75 |
| 5 | −0.61 | −0.72 | ≈ −0.90 | ≈ −0.75 |
| 10 | −0.48 | −0.67 | ≈ −0.65 | ≈ −0.40 |
| 15 | −0.33 | −0.53 | ≈ −0.35 | ≈ −0.20 |
| 20 | −0.26 | −0.38 | ≈ −0.30 | – |

(All values are October of year 15.)

- **Active layer:** 2.0 m (years 1–5), 2.25 m by year 15 with warming. This is consistent with the permafrost table at about 1.9 m in ICT and in the paper (0 °C at about 2.3 m).
- **Below 10 m** the model lies between the paper's numerical and measured curves. **At 3–5 m** it is 0.2–0.4 °C warmer than the paper's curve.
- **Initial-field test.** Starting from a 60-year spin-up without warming (`--spinup 60`) instead of ICT fits better at 10 m but worse at 15–20 m. So ICT is kept. All differences are ≤ 0.4 °C, the same order as the paper's own model-versus-measurement error (0.1–0.2 °C, up to about 0.35 °C at 10 m).
- **Timing.** In this clock the paper's validation day (20 October, surface about +3 °C) is about t = N·8760 − 345 h, two weeks before the month end that `03` reads.
- **Conclusion:** the thermal inputs (properties, Table 1 functions, flux, initial profile) reproduce the paper's natural-ground field to about 0.4 °C.

![1-D temperature profiles](img/column_1d.png)

**Embankment far from the abutment:**
- The permafrost table rises only to about +0.25–0.5 m above the original ground by years 3–5, then drops back to −0.25 m by year 15.
- The fill thaws almost to its base every summer, because the embankment-top mean is +1.5 °C.

### 9.5 Far-field settlement in 1-D, and what the switches do

Far from the abutment (October; total / foundation / embankment, cm):

| Thaw-strain rule | Year 3 | Year 15 |
|---|---|---|
| A: paper (`THAW_REF=0, IRREV=0`) | 11.4 / 4.1 / 7.3 | 12.2 / 4.3 / 7.9 |
| B: irreversible (`IRREV=1`) | 18.8 / 10.8 / 7.9 | 19.0 / 11.0 / 7.9 |
| C: irreversible, ice-poor fill (`IRREV=1, FILL_THAW=0`) | 11.8 / 10.8 / 0.9 | 12.0 / 11.0 / 0.9 |
| D: thawed-at-start strain-free (`THAW_REF=1`) | −3.6 / −3.9 / 0.3 (heave) | −2.8 / −3.7 / 0.9 |
| **Paper, Fig. 14, d = 20 m** | **17.5 / 15.5 / 2.0** | **25.5 / 23.5 / 2.0** |

What this shows:
- **Load compression of the ground is about 4 cm** in every variant. Thaw strain makes the rest.
- **D (heave) contradicts the paper**, so the paper must count thaw strain from the frozen state (A, B or C).
- **The paper's split** (foundation ≫ embankment, body ≈ 2 cm) is closest to **C**. Its year-3 total is closest to **B**.
- **No variant produces the paper's growth from year 3 to 15 at d = 20 m** (+8 cm), because in 1-D the ground far from the abutment does not degrade. In the paper that growth comes from the abutment's heat reaching 20 m (they describe a thawed interlayer "within 20 m") and from 3-D and side-slope effects.

### 9.6 Property and load calculations

- `tools/check_enthalpy.py`: ENTH and latent-heat check (Section 9.3).
- `tools/settlement_check.py` prints every generated E, ν, c, φ, α, σ_Y and creep C₁ table, the creep magnitude, and 1-D thaw settlement against thaw depth. Every 1 m of sub-clay thawed is about **12 cm**; a thaw depth of 6 m under the abutment is about 49 cm of thaw strain alone.
- `tools/traffic_load.py` gives Section 7's numbers.

### 9.7 Not verified without ANSYS: watch these on the first run

| Item | If ANSYS rejects it |
|---|---|
| `TB,EDP` with `TBTEMP` together with `TB,CREEP` (EDP creep) | set `USE_CREEP = 0` (creep is < 1 cm) |
| `MP,THSY` with `MP,REFT` | replace with secant `ALPY` = THSY(T)/(T − T_ref) |
| `MP,DENS` inside `/SOLU` (staged self-weight) | do LS1 in a separate run, or use `EKILL`/`EALIVE` for the fill |
| `*GET,…,NODE,n,BF,TEMP` (only with `IRREV = 1`) | use `*GET,…,NODE,n,NTEMP` |
| `NSEL,R,TEMP` in POST1 (thaw depth) | use `*VGET` of TEMP and loop |

**Look at the log for:**
- the `CNCHECK,SUMMARY` table before LS1: pairs 9 and 10 (bonded) must show contact **closed**;
- convergence in the months when the thaw front is deepest (August to October).

The first real run (ANSYS 18.1) accepted EDP, creep, THSY/REFT and `MP,DENS` in `/SOLU` (the last only with a warning). It stopped on rigid-body motion of the abutment, which is fixed as described in Section 12.1.

---

## 10. What to expect, and the model's limits

### Expected pattern
- Settlement **rises steeply toward the abutment**. The foundation share is **> 80 %**, and the embankment body is a few cm, more next to the abutment.
- Most thaw settlement near the abutment develops in the first years, then grows more slowly.
- The most useful single check is `thaw_depth.csv` near the abutment against the paper (about 3 m in year 3 and 6 m in year 15, April). If it matches, settlement ≈ 4 cm + 12 cm per metre of sub-clay thawed.

### Why the numbers will not match Fig. 14 exactly

1. **Longitudinal 2-D section.** It cannot include the cold side slopes (mean −0.18 °C) or the protection cone, both of which the paper's Table 1 has. So the far-field embankment runs warmer here and the fill thaws more (Section 9.4).
2. The paper does not state its thaw-strain reference, reversibility, creep units, the exponent m, or its model length. The switches cover these; Section 9.5 shows their effect.
3. The paper itself reports 75 cm calculated against about 100 cm measured, and attributes the difference to unknown ice content. Aim to match the **trend and magnitude**, not every point.

### Calibration order (if you want to get closer to the paper)

1. **Thermal first:** get the thaw depth near the abutment right (Fig. 10–12 descriptions). The knobs are the surface T₀ values and the latent heat. Do not tune the mechanics until the thaw depth is right.
2. **Thaw-strain rule:** choose between A, B and C using the embankment-body settlement (paper: about 2 cm far, about 11–12 cm at the abutment).
3. **Stiffness (`M_E`) and load case:** adjust the far-field magnitude.

---

## 11. Parameter reference

**Thermal macro (01)**

| Parameter | Default | Meaning |
|---|---|---|
| `NYEAR` | 15 | years simulated |
| `XL`, `YB` | −50, −30 | left end and bottom of the model (m) |
| `QGEO` | 72 | geothermal flux (J/m²·h) |
| `EMB_MEAN`, `EMB_AMP` | 1.5, 12.5 | embankment-top T₀ and A |
| `FIX_M4` | 1 | restore the mudstone latent peak |
| `USE_CR` (`CR_TOP`, `CR_BOT`, `CR_KS`, `CR_KW`) | 0 | crushed-rock cooling layer |

**Mechanical macro (02)**

| Parameter | Default | Meaning |
|---|---|---|
| `MSTEP` | 1 | months per structural step (1, 2, 3 or 6 keep April and October available) |
| `THAW_REF`, `IRREV`, `FILL_THAW` | 0, 0, 1 | thaw-strain rule (Section 5.2) |
| `THS_1`, `THS_2` | 0, 0.1 | thaw-strain ramp (°C) |
| `VERT_ONLY` | 1 | thaw strain vertical only |
| `USE_CREEP`, `SIG_U`, `TIM_U`, `T_CLMP` | 1, 1000, 1, 0.1 | creep and its units |
| `M_E`, `C_MIN`, `PHI_MX` | 1, 5 kPa, 50° | Eq. 6 exponent and strength safeguards |
| `LOAD_CASE`, `Q_PAVE`, `Q_VEH`, `LOAD_CAP` | 1, –, –, 1 | surcharge (Section 7) |
| `CR_E` | 100 MPa | crushed-rock modulus |
| `ABUT_BOND` | 2 | footing and front: 2 coupled node to node (default), 1 bonded contact, 0 standard (Section 12.2) |
| `FKN_M` | 1.0 | contact normal stiffness factor in the structural run |
| `SEAT_FIX` | 0 | UX = 0 at the bridge seat (superstructure restraint) |
| `STAB` | 1 | energy stabilisation from LS3 on |
| `BEHIND` | 0 | backwall interface: 0 standard (gap), 1 no separation, 2 coupled |
| `USE_PLAS` | 1 | 0 = elastic soil (diagnosis) |
| `MODE_CHECK` | 0 | 1 = quick support check: elastic, all coupled, LS1–LS3 only |

**Extraction macro (03)**

| Parameter | Default | Meaning |
|---|---|---|
| `XA` | −6.2 | abutment back face (d = 0) |
| `ND` | 21 | d = 0 … 20 m |
| `XNAT` | 10 | natural-ground column |
| `YTOL` | 0.26 | column half-width (m) |

---

## 12. How to run, and what to do when it fails

1. Put the three macros in one empty folder. In MAPDL: `/INPUT,01_thermal_transient,mac`. The run is 180 monthly load steps, so it is the long one.
2. Check `POST0*.png` and `thaw_depth.csv` from step 3 part A. (You can run `03` part A alone by stopping it after `FINISH` of part A.)
3. `/INPUT,02_mechanical_settlement,mac` (183 load steps), then `/INPUT,03_extract_results,mac`.
4. Run `python3 tools/compare_fig14.py .`, then `python3 tools/column_1d.py` for the far-field reference.

All the Python checks together (no ANSYS needed):

```
python3 tools/check_geometry.py      # geometry, BC lines, contact, mesh
python3 tools/lint_apdl.py ansys/*.mac
python3 tools/check_enthalpy.py      # ENTH and latent heat
python3 tools/settlement_check.py    # property tables, creep, 1-D thaw settlement
python3 tools/traffic_load.py        # Section 7 numbers
python3 tools/column_1d.py --plot    # 1-D thermal + far-field settlement (~20 s)
python3 tools/compare_fig14.py --demo
```

| Symptom | Fix |
|---|---|
| A month fails to converge | raise the `NSUBST` maximum (e.g. 400), or `MSTEP = 2` |
| Rigid-body error / huge UX, UY at the abutment | see Section 12.1; keep `ABUT_BOND = 1`; optionally `SEAT_FIX = 1` |
| A step yields and will not converge (large plastic zone) | `STAB = 1` (energy stabilisation, 1e-4), check that the stabilisation energy stays small in the output |
| Large settlement already in LS3 | expected in paper mode: the thawed active layer and fill take their thaw strain at LS3 |
| "Sparse solver … out-of-core" warning | only speed. Give MAPDL more memory at launch (e.g. `-m 4000 -db 1024`, or custom memory in the Product Launcher) |
| After an error, "X is not a recognized BEGIN command" | not a separate problem: after a failed SOLVE batch ANSYS leaves `/SOLU` and ignores the remaining solution commands. Fix the first error |
| Settlement jumps every summer and recovers | that is `IRREV = 0`; use `IRREV = 1` for irreversible thaw |

### 12.1 The first run (ANSYS 18.1): what the log meant and what was changed

| Log message | Meaning | Change |
|---|---|---|
| CONTA172 / TARGE169 "no companion element type for ETCHG,TTS" | expected: ETCHG only converts PLANE55 → PLANE182 | none needed; the macro sets the contact KEYOPTs itself |
| "Changing material properties (MP) between load steps is non-standard" | the staged self-weight (`MP,DENS` in LS1 → LS2) was accepted | none |
| Sparse solver "out-of-core" | not enough memory for in-core; slower only | give more memory at launch |
| "N small equation solver pivot terms … had to be constrained", "extremely large pivot ratio" | a part of the model had no stiffness: the **concrete**, held only by contact | see below |
| "UY/UX … greater than 1 000 000 … rigid body motion", "held together only by contact" | in the construction step the concrete got its weight and fell through the soil | see below |
| RATE, CUTCONTR, LDREAD, TIME, NSUBST, SOLVE "not a recognized BEGIN command" | consequence of the failed SOLVE (ANSYS left `/SOLU`) | disappears with the fix |

**Why the concrete had no stiffness:**
1. Concrete and soil have separate nodes. Their surfaces coincide exactly, and standard contact treats "just touching" as **open**, so there is no stiffness.
2. The contact real constants come from the thermal database with **FKN = 0.01**, 100× softer than the default. Even when closed, the 14 m high abutment would sink about 0.2 m into the soil by penetration alone.

**Changes in `02`:**
- `RMODIF` sets FKN = `FKN_M` = 1 on real sets 8–10.
- The footing (real 10) and front (real 9) pairs move to a new element type 4 (CONTA172, `KEYOPT(12)=5`, **bonded always**, augmented Lagrange). The foundation is cast in the ground, so it cannot float or slide away.
- The **behind pair stays standard**, so the embankment can still separate from the backwall, as in the paper.
- `CNCHECK,SUMMARY` prints the contact status before the first solve.
- Loads (LS2) and the first temperature field (LS3) are applied in separate steps; `PRED,ON` is added.
- Optional: `SEAT_FIX = 1` (horizontal restraint from the bridge superstructure at the seat) and `STAB = 1` (energy stabilisation).
- `03` now takes the end of construction from LS3 and counts 3 + months load steps.

### 12.2 Second run: same error. What changed and how to run now

The bonded-contact fix still depends on the contact algorithm detecting the coincident surfaces. The mechanical macro now removes that dependency and reports on itself.

**Changes in `02`:**
- **`ABUT_BOND = 2` (default): footing and front interfaces are coupled node to node** (`CPINTF,UX` / `CPINTF,UY`), and their contact elements are deleted. Every concrete node on those faces has a coincident soil node (same line length and division; checked in Python). So the abutment is joined to the ground exactly, with no contact detection involved.
- **`BEHIND`:** 0 = standard contact (gap allowed, paper), 1 = no separation, 2 = coupled.
- **Supports are taken from the mesh extents** (`*GET … MNLOC/MXLOC`), not from the `XL`/`YB` parameters. The **ICT** profile is redefined in `02`, and `NYEAR` is checked. A missing parameter in `THERMAL.db` can no longer silently move a support.
- **Nodes that belong to no solid element** are counted and fixed (there should be 0).
- A **MODEL CHECK** block prints before solving (search the output for `MODEL CHECK`):
  - model extent;
  - number of fixed nodes on the left, right and bottom;
  - number of nodes, solid elements and free nodes;
  - number of standard, bonded and target contact elements;
  - `CPLIST` of the coupled pairs.
- **`/COM` markers** (`===== LS1 …`, `===== LS2 …`, `===== LS3 …`) show which load step was solving when an error appears.
- **Thaw step LS3:** energy stabilisation (`STAB = 1`, 1e-4) and up to 2000 substeps.
- **`USE_PLAS = 0`** gives elastic soil, to separate plasticity problems from support problems.
- **`MODE_CHECK = 1`:** elastic, no creep, every interface coupled, only LS1–LS3.

**How to run now:**
1. Set `MODE_CHECK = 1` in `02` and run it (a few minutes). It must finish LS1–LS3.
   - If it does not, the problem is outside contact and plasticity. Send the `MODEL CHECK` lines and the node number from the first error.
2. Set `MODE_CHECK = 0` and run the full analysis.
   - If it fails, note which `=====` marker was last printed. LS3 means the thaw step: try `BEHIND = 1`, then `USE_PLAS = 0`, to see which part is unstable.

---

## 13. APDL commands used (glossary)

| Command | What it does here |
|---|---|
| `ET`, `KEYOPT` | element types (PLANE55/182, CONTA172, TARGE169) and their options |
| `MPTEMP`, `MPDATA`, `MP` | temperature tables and constant material properties (C, KXX, ENTH, DENS, EX, PRXY, THSY, REFT, MU) |
| `TB,EDP` / `TBTEMP` / `TBDATA` | Drucker–Prager tables per temperature |
| `TB,CREEP,,,,2` | implicit time-hardening creep |
| `K`, `L`, `AL`, `LESIZE`, `MSHKEY,1`, `AMESH` | geometry and mapped mesh |
| `ESEL,…,CENT` + `MPCHG` | assign material by element centroid |
| `R`, `RMORE`, `ESURF` | contact real constants and generation of contact elements on surfaces |
| `*DIM,…,TABLE` + `%table%` | time-varying surface temperatures (function tables) |
| `D,…,TEMP`, `SF,…,HFLUX`, `IC` | thermal BCs and initial temperature |
| `TIME`, `DELTIM`, `AUTOTS`, `TINTP`, `OUTRES` | transient controls and output frequency |
| `ETCHG,TTS` | thermal-to-structural element change |
| `LDREAD,TEMP` | read temperatures from THERMAL.rth as body loads |
| `BF,…,TEMP` | set nodal body temperature directly (initial profile, irreversible mode) |
| `ACEL`, `SF,…,PRES`, `D,…,UX/UY` | gravity, surcharge, supports |
| `RATE,ON`, `CUTCONTROL,CRPLIMIT` | creep on and the creep-ratio time-step limit |
| `SET`, `*GET`, `UY()`, `TEMP()`, `NODE()` | read results in POST1 |
| `*CFOPEN`/`*VWRITE`/`*CFCLOS` | write CSV files |
| `LCDEF`/`LCOPER,SUB` | subtract the LS1 displacement for contour plots |

---

## 14. Sources

**The paper**
- K. Chen et al. (2024), *Engineering Failure Analysis* 163, 108476: §3.2.2 (Eqs. 4–9), §3.3 (Eq. 10, Tables 1–3, 0.02 W/m², 60.1 kPa), §4 (Figs. 9–14).

**Traffic and pavement load**
- Equivalent characteristic road traffic loads for embankments (15/20 kPa; 31 kPa + 9 kPa model): [Ground Engineering technical paper](https://www.geplus.co.uk/technical-paper/technical-paper-road-traffic-loads-for-geotechnical-analyses-of-embankments-04-06-2020/), [ResearchGate table](https://www.researchgate.net/figure/Equivalent-characteristic-road-traffic-loads-in-kPa-for-plane-strain-conditions_tbl2_346221617), [paper PDF](https://www.researchgate.net/publication/346221617_Road_traffic_loads_for_geotechnical_analyses_of_embankments).
- FHWA/AASHTO traffic surcharge 250 psf (12 kPa): [NYSDOT GDP-11b](https://www.dot.ny.gov/divisions/engineering/technical-services/technical-services-repository/GDP-11b.pdf), [FHWA-HRT-13-046, ch. 6](https://www.fhwa.dot.gov/publications/research/infrastructure/structures/bridge/13046/006.cfm).
- Vehicle load 10 kPa uniform; pavement as 1 m of fill; 30–60 kPa at the expressway subgrade top: [Yang et al. 2021, Adv. Civ. Eng.](https://onlinelibrary.wiley.com/doi/10.1155/2021/9949720).
- Measured dynamic stress 20–50 kPa at the subgrade top, 10–20 kPa at 3.1 m: [An et al. 2018, Shock and Vibration](https://onlinelibrary.wiley.com/doi/10.1155/2018/1956906), [Frontiers 2024, influence depth](https://www.frontiersin.org/journals/built-environment/articles/10.3389/fbuil.2024.1497868/full).
- BZZ-100 (100 kN, 0.70 MPa, d = 21.30 cm, 31.95 cm): JTG D50-2017 ([code page](https://www.codeofchina.com/standard/JTGD50-2017.html)); parameters confirmed in [Baidu Baike: 标准轴载](https://baike.baidu.com/item/%E6%A0%87%E5%87%86%E8%BD%B4%E8%BD%BD/5280650).
- Equivalent soil column for vehicle load in subgrade design: JTG D30-2015 ([code page](https://codeofchina.com/standard/JTGD30-2015.html)); 550 kN heavy vehicle of JTG D60 / JTG B01-2014 ([ch. 7](https://www.erbcc.net/reference/ha-gt/JTGB01-2014/07.html)).
- Eurocode EN 1991-2 Load Model 1 (tandem 300 kN, UDL 9 kN/m²): [JRC bridge design examples](https://eurocodes.jrc.ec.europa.eu/sites/default/files/2022-06/Bridge_Design-Eurocodes-Worked_examples-main_only.pdf).

**Deformation mechanisms on the Qinghai–Tibet Highway**
- [Lateral deformation of expressway embankment on the QTP, IJPE 2023](https://www.tandfonline.com/doi/abs/10.1080/10298436.2023.2171036); [thaw consolidation of permafrost under roadway embankment](https://www.sciencedirect.com/science/article/abs/pii/S0165232X12000882).

Values marked "calculated" in Section 7.2 come from `tools/traffic_load.py`. The JTG vehicle dimensions (L = 12.8 m, b = 1.8 m, m = 1.3 m, d = 0.6 m) are the standard values of those codes; check them against your copy of JTG D30/D60 before citing them in a thesis.
