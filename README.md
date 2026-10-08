# Cargo Aircraft Mission Simulation (MATLAB) Design 27

A time-domain mission simulation for a small electric cargo aircraft: takeoff roll, climb to a cruise altitude, then banked racetrack laps on a 3S LiPo pack. It sits on top of the steady-state sizing tool `cargo_airfoil_optimizer.m` and answers questions the steady-state tool cannot: does the aircraft actually take off, hold altitude in turns, and finish the laps before the battery runs out.

> **Status:** the equations and code have been checked for internal consistency (see [Verification](#verification)), but several inputs are still placeholders or unverified. Treat results as design-trade guidance, not predictions, until the open items at the bottom are closed with measurements.

---

## Files

| File | Purpose |
|---|---|
| `cargo_flight_sim.m` | The simulation (one script, local functions at the end). Needs MATLAB R2019b+. |
| `polars/<airfoil>_Re....csv` | Airfoil polar, columns `Re, alpha_deg, CL, CD`. One block per Reynolds number. |
| `prop/apce_12x6.csv` | Propeller table, columns `J, CT, CP` (optional, see [Propeller data](#propeller-data)). |
| `thrust/thrust_table.csv` | Optional measured thrust-stand data for `cfg.prop.mode = 'table'`. |
| `airfoil_ranking.csv` | Optimizer output (only used if `cfg.useOptimizerResult = true`). |
| `patch_for_current_code.m` | Edit-by-edit patch for an older copy of the script. Not needed if you use the rewritten `cargo_flight_sim.m`. |

Outputs written by the script: `payload_sweep.csv`, `chord_sweep.csv`.

---

## Quick start

1. Put the polar CSV in `polars/` and set `cfg.ac.polarFile`.
2. Set the wing (`span_m`, `chord_m`), masses, battery, motor, and prop in the CONFIG block.
3. Run `cargo_flight_sim.m`.
4. Read the console in this order: **PRE-FLIGHT CHECK**, **MISSION SUMMARY**, **SELF-CHECK**, then the sweeps.

If the pre-flight check warns that the pack hits the cutoff voltage at takeoff, the run will end at t = 0. Lower `cfg.prop.thrMax` or change the pack before reading anything else.

---

## What the model contains

- **Flight dynamics:** 3-DOF point mass, states `x, y, h, V, gamma, psi` plus charge, energy and distance. RK4 integration with zero-order-hold controls.
- **Phases:** (1) takeoff roll with rolling friction and ground effect, rotation at `Vrot = 1.15 * Vstall`; (2) climb at full (capped) throttle; (3) laps: straight legs plus 180 degree turns at a fixed bank angle.
- **Wing:** rectangular, **fixed span and variable chord**. `S = span * chord` and `AR = span / chord` are derived by `deriveWing()`.
- **Aerodynamics:** pre-stall branch of the airfoil polar, interpolated in Reynolds number, plus `CD0_airframe` and induced drag `CL^2 / (pi e AR)` with a ground-effect factor. `CLmax(3D) = 0.9 * CLmax(airfoil)`.
- **External drag articles:** CFD drag forces at two speeds are converted to a drag area `CdA = F / (0.5 rho V^2)` and interpolated in speed.
- **Propulsion:** motor (`Kv, Rm, I0`) coupled to a propeller (`CT(J)`, `CP(J)`) by a closed-form torque balance, with ESC efficiency. `CT` can come from a measured table.
- **Battery:** open-circuit voltage versus state of charge plus series resistance, so voltage sags under load. The mission ends at the reserve SOC or the loaded cell-voltage cutoff.
- **Autopilot:** throttle PI loop on airspeed, altitude-to-flight-path loop, bank-to-heading loop. Cruise speed is `max(Vcruise_ms, VminFactor * Vstall)`.

---

## Configuration reference

| Group | Key fields | Notes |
|---|---|---|
| `cfg.ac` | `span_m`, `chord_m`, `e`, `CD0_airframe`, `CLmax3D_factor`, `CL_clampFrac`, `CL_margin` | `chord_m` is the design variable. `CD0_airframe` is a **placeholder**. |
| `cfg.ac.dragArt` | `V = [7 10]`, `F = [..]` (N) | Total drag of all carried articles at those speeds. |
| `cfg.mass` | `empty_kg`, `dragArticles_kg`, `battery_kg`, `payload_kg`, `wingKgPerM2`, `S_ref` | `wingKgPerM2 = 0` ignores wing weight growth with chord. |
| `cfg.batt` | `cells`, `capacity_Ah`, `R_cell`, `R_wiring`, `reserveSOC`, `cutoffCell_V`, `maxC` | Measure the real DC internal resistance. |
| `cfg.prop` | `mode`, `Kv`, `Rm`, `I0`, `etaESC`, `D_m`, `CT`, `CP`, `CtFile`, `CT0_static`, `thrMax` | `mode` is `'model'` or `'table'` only. |
| `cfg.mis` | `runway_m`, `Vcruise_ms`, `VminFactor`, `nLaps`, `requiredLaps`, `bank_deg`, `legLength_m` | `nLaps = inf` flies to battery reserve. |
| `cfg.sweep` | `payloads`, `payloadNLaps`, `chordStart/End/Step`, `chordNLaps` | Two independent one-variable sweeps. |
| `cfg.check` | `steadyState`, `convergence` | Optional extra checks. |

---

## Input data

### Airfoil polar

CSV with columns `Re, alpha_deg, CL, CD`. Only the branch up to maximum CL is used. Outside the Reynolds range of the file, drag is scaled by `(Re / Re_ref)^ReCDexp` (default -0.5, an assumption). Below the lowest CL in the file, drag is held flat, so digitize the polar down to the CL you cruise at.

### Propeller data

The simulation needs `CT(J)` and `CP(J)`. For the APC 12x6 thin-electric prop, `prop/apce_12x6.csv` is built from the UIUC Propeller Database (Volume 4): the static test supplies the `J = 0` row, and the wind-tunnel run at about 6000 RPM supplies `J = 0.33` to `0.63`.

- Use the **static file** for takeoff and the **highest-RPM wind-tunnel file** for flight. Measured wind-tunnel data does not exist below `J` of about 0.33, so climb thrust is interpolated.
- If the CSV has no `J = 0` row, the script adds one from `cfg.prop.CT0_static` instead of extrapolating (linear extrapolation would overstate static thrust by about 45%).
- `CP` is still a quadratic in `cfg.prop.CP` and drives the RPM solve. Fit values for this CSV: `CP = [0.03366 0.04411 -0.13383]`, `CT = [0.10552 -0.04886 -0.19250]`.
- Keep `CT` and `CP` from the same propeller. For other sizes use APC performance files or another measured source.

### Motor

`Kv`, `Rm` and `I0` from the datasheet. Confirm the voltage `I0` was measured at.

### Drag articles (CFD, isolated bodies at 0 degrees angle of attack)

| Article | Min. mass | Drag at 7 m/s | Drag at 10 m/s | CdA (m^2) |
|---|---|---|---|---|
| Waffles | 226 g | 0.8385 N | 1.6740 N | 0.0279 |
| Dumpster | 285 g | 0.4717 N | 0.9561 N | 0.0157 |
| Flying pig | 340 g | 0.2675 N | 0.5429 N | 0.0089 |
| Banana | 453 g | 0.1496 N | 0.3020 N | 0.0050 |
| Chicken | 680 g | 0.0241 N | 0.0490 N | 0.0008 |

Exactly three are carried, all external. For Pig + Banana + Chicken: 1.473 kg, `F = [0.4412 0.8939]` N. In the simulations run so far this combination gave the best range. Mounted on the aircraft the drag may be higher than the isolated CFD value (wing wake, prop slipstream, interference), so run a sensitivity case.

---

## Outputs

- **Console:** pre-flight check, mission summary (status, ground roll, laps, distance, energy, peak current, minimum cell voltage, time near stall), self-check, steady-state answer key.
- **Figures:** flight (ground track, altitude, airspeed, CL), battery and propulsion, power versus distance, Reynolds number, power versus time, payload sweep, chord sweep, steady-state power.
- **Mission status values:**
  - `COMPLETE: commanded laps flown`
  - `END: battery reserve SOC reached`
  - `END: loaded cell voltage below cutoff`
  - `FAIL: runway exceeded before Vrot`
  - `FAIL: ground impact`

---

## Sweeps

- **Payload sweep:** each payload is flown at the configured chord for `payloadNLaps` laps. The plot shows charge used and ground roll.
- **Chord sweep:** each chord is flown at the configured payload for `chordNLaps` laps (`inf` = to battery reserve). A chord is marked acceptable only if the mission passes, turn CL stays below `CL_margin * CLmax`, and `Vrot` is below the cruise command.
- The two sweeps are independent. A chord-by-payload grid is not run.
- Failed runs report NaN distance. The distance of a crashed run is only how far it got and must not be compared with a completed run.

---

## Verification

Built-in checks:

- **`preflightCheck`:** static thrust, current, and pack voltage at the throttle cap. Warns if the run will end at t = 0.
- **`selfCheck`:** integrated power and current match the energy and charge state; thrust equals drag and load factor is 1.00 in steady cruise; load factor is `1/cos(bank)` in turns; altitude is held; lift is not pinned at the clamp.
- **`steadyStateCheck`:** steady-level-flight "answer key": stall speed, L/D, power required, electrical power, best endurance and range speeds, maximum speed, takeoff estimate, and the equivalent parabolic `CD0`.
- **`convergenceCheck`:** reruns at `dt = 0.04, 0.02, 0.01`. Differences below about 1% between 0.02 and 0.01 mean the step is fine.

Cross-checks outside the script:

- **MathWorks tools** (Aircraft Performance Analyzer, Constraint Analysis, Aircraft Intuitive Design): enter the same weight, area, span, efficiency, equivalent `CD0`, `CLmax`, overall propulsive efficiency and usable capacity, then compare stall speed, power required, endurance and takeoff roll. Align assumptions first, then interpret what is left.
- **Hardware tests**, which are the only true validation: thrust stand for static thrust and current, glide test for `CD0_airframe`, stall test for `CLmax`, DC internal resistance from the charger, and logged cruise power against the simulation.

Agreement between tools shows consistency, not truth, because they share the same inputs.

---

## Known limitations

- Point-mass model: no CG, static margin, tail sizing, control surfaces, roll dynamics, or stall behaviour (CL is clamped at 0.95 of CLmax rather than stalling).
- Rectangular wing only. Span efficiency `e` is constant across aspect ratio. Wing weight does not grow with area unless `wingKgPerM2` is set.
- Drag at a single polar Reynolds number is scaled with an assumed exponent outside the data.
- CFD drag articles are isolated bodies at zero angle of attack with no interference.
- Takeoff uses a fixed ground-roll CL and rolling friction; the runway length is a placeholder.
- Energy for the climb is large: with full-throttle climbing, takeoff and climb used about half the charge of a two-lap mission in earlier runs.

---

## Troubleshooting

| Symptom | Likely cause |
|---|---|
| `END: loaded cell voltage below cutoff` at t = 0 | Static current too high for the pack. Lower `thrMax`, use a bigger or lower-resistance pack. |
| `FAIL: ground impact` soon after the first turn | Wing too small or too heavy for the cruise speed. Raise chord, lower payload, or raise `VminFactor`. |
| `FAIL: runway exceeded before Vrot` | Runway too short for the weight and thrust. Check `runway_m` against the rules. |
| Distance jumps between neighbouring sweep points | One of them crashed. Look at the status column. |
| `cfg.prop.mode must be 'model' or 'table'` | Any other mode string. |
| `Prop CT file not found` | Fix `cfg.prop.CtFile` or set it to `''`. |
| Throttle stuck at 100% in the laps | The speed-controller throttle line was removed or overwritten. |

---

## Open items

1. `CD0_airframe` is a placeholder. It is the largest unknown affecting range.
2. Confirm motor `Kv`, `Rm` and `I0` against the actual datasheet.
3. Measure the pack's DC internal resistance (and confirm the real continuous C rating).
4. Review the polar for the chosen airfoil down to the cruise CL, and confirm its `CLmax`.
5. Set the real runway length from the rules. Check that payload and runway are consistent.
6. Check the name and file pairs in `applyOptimizerResult` before enabling the optimizer.
7. Compare against MathWorks tools and run the hardware tests above.

---

## Data sources

- Propeller coefficients: UIUC Propeller Database, Volume 4 (APC 12x6 thin electric).
- Drag articles: the team's own Ansys Fluent analysis (k-omega SST with gamma transition model, 7 and 10 m/s).
- Airfoil polars: digitized from published polars at the stated Reynolds numbers.
