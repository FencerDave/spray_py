# Wish List, Known Issues & Design Notes

Living document for **spray_py**: what the model deliberately ignores, what we would
like to add, and the known loose ends. Compiled from the 2026-10 design/review sessions.

Units throughout: **m, g, J, s, °C** (temperature differences only, so no Kelvin needed).

---

## 1. Deliberately out of scope (v1)

These are real physics, and a big world of its own. Left for a future contributor
(colleague, academic) on purpose, not by oversight.

- Shrinkage, crust formation, hollow shells, rearrangement of solids (microstructure physics)
- Porosity-dependent conductivity and density (dry material currently keeps the same
  volume; `dry_condK` is a constant)
- Evaporative cooling and vapor transport at the surface (the model is conduction plus a
  latent-heat sink only; real droplets are limited toward the wet-bulb temperature)

---

## 2. Physics improvements

- **Air temperature profile along the dryer.** `sprayType`: Constant (current), linear
  inlet to outlet (rough co-current), dual-linear or quadratic (rough counter-current).
- **Residence time** from drag, particle velocity and swirl instead of a fixed `dwellTime`;
  or run until the droplet reaches e.g. 95% of the air temperature.
- **Convection coefficient `h`** from a Nusselt correlation (Ranz-Marshall) instead of a
  constant 1000 W/m²K. A stationary droplet (Nu = 2) gives roughly 100-200 W/m²K at
  600 µm diameter; the droplets here are not stationary. `h` sets the Biot number, so sweep it.
- **Outlet air energy balance** and humidity (`InletAir.flowRate`, `moisturePct`,
  `heatCapacity` are currently unused).
- **Mixture heat capacity** for the wet state (water + solids), instead of water's Cp
  applied to the whole voxel mass.
- **Interface conductivity at wet/dry boundaries**: harmonic mean `2·k1·k2/(k1+k2)`.
  Conductivity across the boiling region (currently switches only when fully dry).
- **Droplet size distribution** (`dropR_sigmaDist`): model a larger and a smaller droplet
  to check whether the largest are undried or the smallest are burnt.
- **Recipe-based properties**: estimate k, Cp and density from `solvent` / `solids`.
- Thermodynamic constants collection ("Universe" dataclass), if actually needed.

## 3. Numerics

- **Phase-boundary look-ahead with a floor.** Do not let one step carry a voxel across
  boiling to dry with lots of leftover energy (it lands on a ~2800x smaller heat capacity
  and kicks the temperature). Limit `dt` by remaining enthalpy budget / net power **with a
  floor**, or let the step cross the boundary and carry the remainder into the next phase.
  Without a floor, `dt` shrinks toward zero as a voxel nears the boundary (Zeno effect).
  Equivalent to the old "constant heat packet per step" method from the notebook.
- **Check ringing after the whole droplet is dry** (surface overshooting the air
  temperature). Physical invariant: `T <= T_air` always.
- **Semi-implicit surface convection** as a conservative backstop:
  `dQ = hA(T_air - T)·dt / (1 + hA·dt/C)`.
- **Implicit time stepping** (tridiagonal backward Euler) or `scipy.integrate.solve_ivp`
  with a stiff method and events for phase changes, so `dt` is not stability-limited.
- **Grid-convergence study** to document the choice of `dR` (e.g. 10, 5, 2.5 µm).
- Raise `dt_max` (currently 1 µs) when the wet-phase limit allows; at `dR = 5 µm` the cap,
  not stability, was the binding limit.
- Dynamic `dt` should clip to land exactly on log times and the end time.

## 4. Analysis & workflow

- History recorder with preallocated arrays, converted to a DataFrame at the end
  (separate from the solver).
- "Values of note" reporting: mass-averaged T and moisture, surface and center T, % mass
  under the target and over the damage temperature.
- **Pass/fail map**: (A) center above the target temperature, (B) surface above the damage
  temperature. A and not B = success; A and B = overheated, lower the air temperature;
  not A and not B = underdried, raise it; not A and B = failure (droplet size and dryer mismatch).
- Plots: T(r) lines per time with bold every Nth; surface/center/average vs time;
  pass/fail colouring.
- Precompute an expected time-to-dry from an energy balance, to size runs
  instead of guessing `dt` and `N`.
- Parameter sweeps and optimization (inlet temperature, flow, solids %).
- Minimum viable **Streamlit app** front end.

---

## 5. Known issues / cleanup (non-blocking)

| Item | Where | Note |
|---|---|---|
| Convection term in stability limit | `_stable_dtime`, `sigG_vox[-1]` | Uses `k_upper[-1]*a_upper[-1]` (~4e-9); should be `h*a_upper[-1]` (~2e-4). Negligible numerically, but wrong. `h` is not passed in. |
| `geoAlpha = 3` double-counts geometry | `_stable_dtime` | `C/ΣG` already includes the geometry (center voxel gives `dR²/(3α)`). Use `geoAlpha = 1` with `margin = 1.0` (dt/2). |
| Wrong unit label in `dt_min` / `dt_max` warnings | `_stable_dtime` | `round(dt*10**6)` is microseconds, labelled "nanosec". |
| `dt_check` default of `0` as a sentinel | `_stable_dtime` | Use `None` for "don't check". |
| Recomputing the limit every step when `dt` is fixed | `iterate` | Compute once and after phase changes. |
| Docstring says `dt_lim / geometry * (1 + margin)` | `_stable_dtime` | Code divides by `(1 + margin)`. |
| Two sources of truth for `dR` | `SprayDryer` vs `SprayDroplet(dR=...)` | Keep one. |
| Dead attributes | `dT_upper`, `dT_lower`, `isDry`, commented loops | Delete. |
| Builtin `max()` / `min()` / `sum()` on arrays | `iterate`, loop | Use `.max()`, `.min()`, `.sum()` (the old `sum(a, b)` bug came from this). |
| Volume-weighted solids average | `solidsPct_avg` | Use `m_solid.sum() / m_voxel.sum()`. |
| Conductivity unit labels | dataclass comments | Thermal conductivity is W/(m·K), not "Watts*m/m2 K". |
| Comments | `inletTemp  # *C--> K` | Temperatures are °C throughout. |
| camelCase attribute names | dataclasses | Decide a PEP 8 migration (snake_case). |
| Module-level script code in the library file | `UnitOps.py` | Import side effects block pytest. Split library from driver (first task of next session). |
| `0 * np.inf` gives NaN | any boolean-mask math with `inf` | Use `np.where`. |

---

## 6. Test list (to build, with `pytest`)

1. `_stable_dtime` returns a finite positive number in every phase (no NaN)
2. Equals an independent `min(C/ΣG)` calculation, including the `h·A` surface term
3. A deliberately huge `dt` warns; a healthy `dt` does not
4. `Age` equals the sum of the `dt` values used (dynamic-dt bookkeeping)
5. Mesh: outer edge equals `dropR`; Σ voxel volume equals 4/3πR³ for several `dR`, including non-dividing ones
6. Energy ledger stays at 100% (tight tolerance) through heating and boiling
7. Invariants: `T <= T_air`, `T >= T_init`, solids % non-decreasing and within [init, 1]
8. Uniform T with air at the same T stays uniform
9. Low-Bi lumped-capacitance comparison vs `exp(-t/τ)`, `τ = ρ·Cp·R/(3h)`
10. `dt` convergence (halve `dt`, small change); dynamic `dt` ≈ small fixed `dt`
11. One golden-master regression run (a few key numbers)

## 7. Debugging practices to adopt

- Model returns data; plotting is a separate step.
- Physics checks as debugger: energy ledger and invariants (optional debug mode).
- `np.seterr(all="raise")` to stop at the first NaN or overflow.
- Spyder: conditional breakpoints (Shift+F12), `%debug` after a crash, `#%%` cells.
- `logging` instead of `print`; `warnings.warn` once, not per step.
- Log the first violation with context (`dt`, voxel index, phase, limit) and stop.

---

## 8. Design notes (what we established, and why)

- **Enthalpy is the state variable** (`Qnet_vox`). Temperature and solids % are derived from it
  piecewise: heating slurry, boiling plateau (`Hvap`), heating dry solid. This removes the
  edge cases of overshooting 100 °C or leftover heat on the step where boiling ends.
- **Conservative flux**: compute each interface flux once, add it to one neighbour and
  subtract it from the other (`q[:-1] += flux; q[1:] -= flux`). The first version used
  `np.gradient` per voxel, which leaked energy (stored/in was only 11-39%) and caused the
  "stuck at ~100.8 °C" plateau.
- **Conservation is not stability.** An unstable `dt` still conserved energy but zig-zagged
  every step. Test both.
- **Explicit stability limit** per voxel: `dt < C / ΣG` (this is the no-ringing bound).
  True threshold measured at about 1.46x higher (eigenvalue limit). Center voxel binds.
  Dry material has about 50x the diffusivity of wet, so dry sets the worst case.
  At `dR = 5 µm` (dry): bound 1.33 µs, exact limit 1.94 µs.
- **Biot number** (L = R/3): about 0.2 wet, 0.04 dry at `h = 1000`, R = 300 µm. Whole-droplet
  Bi sets how steep gradients are; the grid Biot `h·dR/k` (~0.002) shows convection barely
  affects `dt`. Cost scales as 1/dR³.
- **Mesh**: voxels identified by lower radius; outer edge must equal `dropR`.

## 9. References (verify editions and section numbers before citing)

- Incropera, DeWitt, Bergman & Lavine, *Fundamentals of Heat and Mass Transfer*, ch. 5
  (transient conduction: lumped capacitance Bi < 0.1; explicit finite-difference stability).
- Çengel & Ghajar, *Heat and Mass Transfer: Fundamentals and Applications* (transient and
  numerical-methods chapters).
- Carslaw & Jaeger, *Conduction of Heat in Solids* (analytic sphere solution with a convective boundary).
- Patankar, *Numerical Heat Transfer and Fluid Flow* (1980) (conservative control-volume formulation).
- Alexiades & Solomon, *Mathematical Modeling of Melting and Freezing Processes*; Crank,
  *Free and Moving Boundary Problems* (enthalpy method for phase change).
- Ranz & Marshall (1952), *Chem. Eng. Prog.* 48 (droplet heat transfer coefficient).
- Masters, *Spray Drying Handbook* (spray-dryer design context).
