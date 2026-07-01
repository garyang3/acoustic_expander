# Acoustic Expander — project notes

Computational-physics study of an **acoustic expander**: a gas oscillator (Helmholtz-style
neck + compression volume + passive volume) driven through a **reed valve**. The work explores
when the oscillation is self-sustaining ("success") vs. decays ("fail") across parameter sweeps,
plus energy-balance and oscillation-envelope analysis. Mathematica `.nb` files cover the reed model.

## Environment

- Python via conda/miniforge: `/Users/a1/miniforge3/bin/python` (use this interpreter, already allowed in settings).
- Core stack: **numpy**, **scipy** (ODE integration), **matplotlib** (PDF/PNG output), **joblib** (parallel sweeps). Occasional `pandas`, `gc` (free memory between large sweeps).
- Work is in **Jupyter notebooks** (`.ipynb`), not a package. There is no test suite, requirements file, or build step.

## Notebook map

- `Acoustic_Expander.ipynb` / `Acoustic_Expander_22.ipynb` — baseline model + the main 2D sweeps (P0×v0, valve size, V0_cv×V0_pass, m_gas×V0_pass) producing heatmaps.
- `time_varying_Pin.ipynb` — ramped inlet pressure `P_in_ramp(t, lam)`; effect of λ and `m_max`.
- `V_passive_sweep_energy_balance.ipynb` — enthalpy-in vs. damping-loss energy balance.
- `zeta_scan.ipynb`, `m_max_scan.ipynb`, `m_max_zeta_heatmap.ipynb` — single-parameter scans.
- `fail_vs_success_analysis.ipynb`, `trajectory*.ipynb`, `*amplitude*` — classify/visualize outcomes.
- `linear_modeling.ipynb`, `scplot.ipynb` — linearized analysis / unrelated superconducting-memory plots.

## Physical model (recurring symbols — don't re-derive each session)

- Gas: `gamma=1.41`, `R=287`, `cp=1050`, `cv=744` (J/kg·K), `T_in=300 K`.
- Geometry/inertia: `A_neck=2.85e-4 m²`, `L_neck=0.15 m`, `m_gas≈0.0015 kg`, `V0_cv` (compression vol.), `V0_pass=3.93e-4 m³` (passive vol.).
- Pressures: `P_in≈6e5`, `P_out≈5.45e5`, `P0≈5.725e5` Pa (initial average).
- Drive/damping: `f≈130–135 Hz`, quality factor `Q`, damping `b_gas` (note `b_gas = omega*m_gas/Q`), `zeta`, reed cap `m_max`.
- Reed valve: `reed_flow(dP, m_max, dP_max)` — mass flow saturating at `m_max`.
- ODE right-hand sides: `make_ode(...)` / `acoustic_expander_ode[_asymmetric](t, y)`; sweeps run `simulate_once(...)` / `_sim_*` helpers in parallel via joblib `_*_task` wrappers.

## Conventions & gotchas

- **Notebooks are large (0.4–2.3 MB).** Don't dump full notebook JSON — extract the cells you need with a small Python script (parse `nb["cells"]`), or use the `notebook-explorer` subagent.
- Plot outputs land in `plots/`, `Oscillation_Plots/`, and per-case dirs (`different_valve_size/`, `v=-16/`, `v=0.1/`, `middle_pressure_different_valve/`). Filenames encode parameters, e.g. `..._P_0=578225.0__v0=10.00__f=130.pdf`.
- Large sweep results are cached as `.npz` (e.g. `run_dP_max_out=1.78.npz`, ~32 MB) — load rather than recompute.
- Constants are redefined per notebook and drift slightly between them (e.g. `V0_cv` differs); trust the values in the notebook you're editing, not another.
