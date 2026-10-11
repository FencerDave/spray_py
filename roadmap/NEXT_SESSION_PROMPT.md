# Next-session prompt: project hygiene + test-driven setup

How to use: start a new Claude Code session in this project and paste everything below the line.
Commit and push first, so the session sees your latest `spraypy/UnitOps.py`.
Companion document: `roadmap/WISHLIST.md` (known issues, test list, design notes, references).

---

# Project: spray_py — spray-dryer droplet model. Goal now: project hygiene + test-driven setup

## Who / how to work with me
I'm David (dhooke), a chemical engineer returning to Python after ~5 years. I want
hand-holding and explanations for tooling (pytest, ruff, coverage, CI) — explain what each
command does and why. **I write the model code myself; you advise and review.** Don't edit
my model code unless I ask. Tooling/config files (pyproject.toml, tests scaffolding) are fine
for you to create when I ask, but explain them. One step at a time. I commit/push at
stopping points; if you work in a git worktree, sync from origin/master before reading code.
Read `roadmap/WISHLIST.md` first: it holds the known-issue list, test list and design notes.

## Project in one paragraph
Models one spray-dried droplet's life: a sphere discretized into radial shells (voxels).
Convective heat from hot air enters the outer shell, conducts inward (explicit Euler,
forward time stepping), and each voxel passes slurry → boiling plateau (latent heat of
vaporization) → dry solid (new k, Cp). **State variable is cumulative enthalpy per voxel
(Qnet_vox)**; temperature and solids % are derived from it by piecewise logic. Units: m, g, J,
s, °C. Repo: github.com/FencerDave/spray_py (branch master). Legacy notebook/script now in
`legacy/`. Only code file: `spraypy/UnitOps.py` (~650 lines): dataclasses `InletAir`,
`SlurryData`, `SprayDryer`; class `SprayDroplet` (`__init__`, `_stable_dtime`, `iterate`);
then ~120 lines of module-level driver/plot code (the loop, matplotlib figures, timing).

## What's already been established (verified by running it)
- Conduction is conservative (one flux per interface, N+1 array math); energy ledger
  (sum(Qnet_vox)/Qnet_drop) = 100.000% in runs.
- Stability: explicit limit per voxel dt < C / ΣG (C = voxel heat capacity J/K, ΣG = sum of
  conductances to neighbours + h·A at the surface). That's the *no-ringing* bound. True
  instability threshold (eigenvalue) is ~1.46x higher. Center voxel is the binding one
  (dR²/(3α)). Boiling voxels: C=inf (skip). Dry α ≈ 50x wet α, so dry limits dt.
  Conservation ≠ stability: an unstable dt still conserved energy but zig-zagged every step.
- Biot number (L=R/3): ~0.2 wet, ~0.04 dry at h=1000, R=300µm. Grid Biot h·dR/k ~0.002 so
  convection barely affects dt. Cost ∝ 1/dR³. dR choice must be justified by a
  grid-convergence study, not by Bi directly.
- `iterate(..., dynamic_dt=True)` works: dt now finite, Age tracks true time, warnings are
  not silent. Run at dR=5µm sat at dt_max=1µs for 100 ms; dR=1µm dt ≈ 0.2–0.6µs.
- NOT yet verified: behavior after the whole droplet is dry (earlier concern: dry outer voxel
  overshooting air temperature = ringing). Physical invariant: T must never exceed T_air.

## Known issues to fix later (not now) — see roadmap/WISHLIST.md section 5
- Line ~245 `sigG_vox[-1] += k_upper[-1]*a_upper[-1]` should be `h*a_upper[-1]` (h not passed in).
- `geoAlpha=3` double-counts geometry (C/ΣG already includes it); use 1, with margin=1.0 (= dt/2).
- Phase-boundary look-ahead with a floor (avoid Zeno shrinking dt / big kick when a voxel
  flips boiling→dry); alternatively "max enthalpy packet per step" like my old notebook.
- Semi-implicit surface convection as backstop; implicit/solve_ivp long-term.
- Harmonic-mean k at wet/dry interfaces; mass-weighted solids average (currently volume);
  mixture Cp for the wet state; Ranz–Marshall h; dt_min warning prints wrong unit label.
- Dead attributes (dT_upper, dT_lower, isDry, dropR/dRadius duplicates), `max()` builtin on arrays.
- Out of scope on purpose: shrinkage/hollow particles/microstructure, evaporative surface
  cooling, porosity effects (recorded in WISHLIST.md "deliberately out of scope").

## Environment facts
Windows 11, Spyder, Anaconda base (Python 3.12.7) has pytest 7.4.4, streamlit 1.37.1, black,
mypy, numpy 2.3.4, matplotlib 3.9.2 — BUT matplotlib crashes under numpy 2.x there; no ruff,
no pytest-cov. Want a clean project environment (venv or conda env) that Spyder can also use.

## Plan (do in this order; confirm each step with me)
1. **Environment + packaging**: pyproject.toml (project metadata, dependencies numpy,
   matplotlib; optional extras `dev` = pytest, pytest-cov, ruff; `app` = streamlit), clean env,
   `pip install -e .[dev]`, `spraypy/__init__.py`, `tests/` folder, update .gitignore. Make
   Spyder use the same interpreter. Explain editable installs.
2. **Remove import side effects**: split `UnitOps.py` into a library (dataclasses +
   SprayDroplet, no plotting/printing/top-level code) and separate driver/plot code
   (e.g. `run_demo.py` with `if __name__ == "__main__":`, plotting in its own module).
   I'll propose names; help me decide the module layout. Keep behavior identical (verify
   with the energy ledger + a saved reference result).
3. **First tests** (I write them, you review; start tiny, with a fast fixture: dR=10µm):
   a. `_stable_dtime` returns a finite positive number in all phases (no NaN)
   b. equals an independent min(C/ΣG) computation, including the h·A surface term
   c. a deliberately huge dt warns; a healthy dt does not (pytest.warns / warnings-as-errors)
   d. Age == sum of dt actually used (dynamic dt bookkeeping)
   e. mesh: outer edge == dropR; Σ voxel volume == 4/3πR³ for several dR incl. non-dividing
   f. energy ledger stays 100% (tight tolerance) through heating AND boiling
   g. invariants: T ≤ T_air, T ≥ T_init, solids % non-decreasing, within [init, 1]
   h. uniform T with air at same T stays uniform (conduction sanity)
   i. low-Bi lumped-capacitance comparison vs exp(-t/τ), τ = ρ·Cp·R/(3h)
   j. dt-convergence (halve dt → small change) and dynamic-dt ≈ fixed-small-dt
   k. one golden-master regression run (store a few key numbers)
   Use pytest basics: fixtures, parametrize, `pytest.approx`, markers (`slow`), and how to
   read failures; show the Spyder/terminal workflow.
4. **Quality tools**: ruff (lint + format, replacing black), maybe pre-commit, type hints
   basics (mypy optional), PEP 8 naming (my attributes are camelCase — decide a migration).
5. **Coverage + CI**: pytest-cov, GitHub Actions workflow (run tests + ruff on push/PR),
   coverage badge in README (Codecov or similar) — I remember doing this years ago but need
   it walked through.
6. **Debugging practices**: logging instead of print, `%debug` / post-mortem in Spyder,
   conditional breakpoints (Shift+F12), `np.seterr(all="raise")`, invariant checks as an
   optional debug mode.
7. **Minimum-viable Streamlit app** (`app/` folder): sliders for air temp, h, droplet radius,
   solids %, dR/dt; run the model through the library (no side effects on import), plot T(r)
   and solids(r) at chosen times plus center/surface vs time; use `st.cache_data`. Goal: prove
   the library imports and runs cleanly from a front end.
8. Then return to the physics TODO list above and the WISHLIST.md.

## First thing to do in this session
Look at the repo (`git log`, `spraypy/UnitOps.py`), check my environment, and walk me through
steps 1–2 before any tests. Ask me questions where a choice is mine.
