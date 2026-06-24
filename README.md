# Ocen IMBH Mock Analysis

This repository is a mock-data companion to `Ocen_IMBH_analysis`. Which contains scripts for fitting orbits of high-velocity stars around an intermediate-mass black hole (IMBH) candidate in ω Centauri. It adds utilties to generate mock 
single-epoch position, proper-motion, acceleration, and radial-velocity observations 
from input orbital parameters and then fits them with the same pipeline used for the real data.

The repository can still be used to fit the real `octo_utils.stars` dictionary by setting
`[mock].enabled = false`.

Orbit fitting uses [`Octofitter_imbh.jl`](https://github.com/vincent-hb/Octofitter_imbh.jl.git), a development fork of [Octofitter](https://github.com/sefffal/Octofitter.jl), sampling with Pigeons (parallel tempering HMC/NUTS).

---

## What Changed from [`Ocen_IMBH_analysis`](https://github.com/vincent-hb/Ocen_IMBH_analysis.git)

- Added mock data generation helpers in `launch_scripts/octo_utils.jl`:
  `make_star`, `mock_data`, `stardata_struct`, and `build_mock_observations`.
- Extended `StarData` for position uncertainties
  (`sigma_ra`, `sigma_dec`) so mock position errors can be controlled from TOML.
- Added a `[mock]` TOML section containing input IMBH parameters, observation
  uncertainties, and orbital elements for each star.
- Updated `octo_orbit_direct_likelihoods.jl` to check `[mock].enabled`.
  When enabled, it builds mock-star observatoins from the mock config instead of using
  the real-data dictionary.
- Added scripts for analysis and forcasting:
  `mock_checks.jl`, `forecast_curvature.jl`, `chain_samples.jl`,
  `summarize.jl`, and `batch_posterior_plots.jl`.

---

## Repository Structure

```text
Ocen_IMBH_Mock_Analysis/
|-- configs/
|   |-- mock_MAP.toml
|   |-- mock_max_logpost_closest_mass.toml
|   |-- mock_rand_1.toml
|   |-- mock_rand_2.toml
|   `-- mock_rand_3.toml
|-- launch_scripts/
|   |-- octo_orbit_direct_likelihoods.jl   # Main fitting script
|   |-- octo_utils.jl                      # Utilities 
|   |-- parse_config.jl                    # TOML loading, prior/data helpers
|   |-- plot_chain.jl                      # Posterior and orbit plots
|   |-- submit_job.jl                      # Slurm job generator/submitter
|   |-- mock_checks.jl                     # Compare mock orbits to real positions (debugging)
|   |-- forecast_curvature.jl              # Future curvature detectability
|   |-- chain_samples.jl                   # Chain extractor
|   |-- summarize.jl                       # Summarize many runs 
|   |-- batch_posterior_plots.jl           # Batch plotting script
`-- mock_results/  *full directory is too large for GitHub    
    |-- run_outputs/                       # Chains, summaries, plots
    `-- logs/                              # Slurm stdout/stderr and generated jobs
```

---

## Prerequisites

- Julia 1.10+ or the cluster module specified in the config
  (currently examples use `julia/1.11.3`).
- Octofitter IMBH environment with `../Octofitter_imbh.jl`:

`submit_job.jl` runs Julia with the `[paths].project` environment and
adds some plotting and sampling dependencies before launching the fit.
The main fitting script also ensures `OctofitterRadialVelocity` is developed
from `Octofitter_imbh.jl`.

*** Note that not all dependencies are automatically loaded. Sorry! (At least KernalDensity is missing)

---

## Quick Start

Run a mock fit locally or in an interactive session:

```bash
julia --project=../octoIMBH_env 
    launch_scripts/octo_orbit_direct_likelihoods.jl \
    configs/mock_MAP.toml
```

From inside `launch_scripts/` on the cluster:

```bash
julia --project=../../octoIMBH_env octo_orbit_direct_likelihoods.jl ../configs/mock_MAP.toml
```

Submit through Slurm:

```bash
cd launch_scripts/
julia submit_job.jl ../configs/mock_MAP.toml --dry-run
julia submit_job.jl ../configs/mock_MAP.toml
```

The dry run writes a timestamped Slurm script to `mock_results/logs/` without
submitting it.

---

## Configuration Files

All run parameters live in TOML files under `configs/`. 

### `[meta]`

Human-readable labels written into the run summary.

```toml
[meta]
system_name = "Mock_Omega_Cen"
description = "Mock dataset 5-star direct likelihood fit"
```

### `[stars]`

Which stars to include. Existing labels are `A` through `G`; most mock configs
use `A`, `C`, `D`, `E`, and `F`.

```toml
[stars]
selected = ["A", "C", "D", "E", "F"]
```

### `[epoch]`

Reference epoch at which position, proper motion, acceleration, and RV are
evaluated. The decimal year is converted to MJD internally.

```toml
[epoch]
year = 2010.0
```

### `[priors.system]`

System-level priors shared across all stars.

| Parameter | Meaning | Units |
|-----------|---------|-------|
| `plx` | Parallax | mas |
| `M` | IMBH mass | solar masses |
| `offsetx` | IMBH RA* offset from assumed cluster center | mas |
| `offsety` | IMBH Dec offset from assumed cluster center | mas |

```toml
[priors.system]
plx     = "truncated(Normal(0.19, 0.004), lower=0)"
M       = "Uniform(100, 100000)"
offsetx = "Uniform(-3000, 3000)"
offsety = "Uniform(-3000, 3000)"
```

### `[priors.companion_defaults]`

Default orbital priors applied to every selected star unless overridden.

| Parameter | Meaning |
|-----------|---------|
| `P` | Period [yr] |
| `e` | Eccentricity |
| `i` | Inclination; use `"Sine()"` |
| `omega` | Argument of periastron; use `"UniformCircular()"` |
| `Omega` | Longitude of ascending node; use `"UniformCircular()"` |
| `theta` | Mean anomaly at the reference epoch; use `"UniformCircular()"` |

Optional per-star overrides go in `[priors.overrides.<STAR>]`.

```toml
[priors.overrides.A]
P = "Uniform(10, 5000)"
```

Supported prior strings:

```text
"Uniform(lo, hi)"
"Normal(mu, sigma)"
"truncated(Normal(mu, sigma), lower=L)"
"truncated(Normal(mu, sigma), lower=L, upper=U)"
"Sine()"
"UniformCircular()"
```

Underscores in numbers are allowed, for example `"Uniform(10, 2_000_000)"`.

### `[data.defaults]` and `[data.overrides.<STAR>]`

Controls which observations enter the likelihood. Defaults apply to all stars;
per-star overrides disable or enable individual data types.

```toml
[data.defaults]
position        = true
proper_motion   = true
acceleration    = false
radial_velocity = true
escape_velocity = false
position_oneil  = true

[data.overrides.A]
radial_velocity = false
```

Notes:

- `position_oneil = true` wraps position observations in
  `ObsPriorAstromONeil2019` prior.
- Radial velocity is only included when the generated or real `StarData` has RV
  data; current configs disable it for A, C, and D.
- The z prior is controlled separately by `[data.z_prior]`.
- `escape_velocity` adds the Haberle-style piecewise escape-velocity
  constraint when enabled.

### `[data.z_prior]`

Adds a line-of-sight prior on each star's z coordinate.

```toml
[data.z_prior]
sigma_z_au = 4558.0
```

The prior is `Normal(0, sigma_z_au)` in AU. Some configs use `4558 AU`, based on
the one-dimensional positional standard deviation for the fast stars; this could 
also be taken as the Omega Cen core radius (`845000 AU`).

### `[mock]`

Turns mock-data generation on or off and defines the input parameters.

```toml
[mock]
enabled = true

M_IMBH  = 77306.400872693033
plx     = 0.187909631720
offsetx = -122.205858303562
offsety = 779.597275235741

sigma_ra      = 0.5
sigma_dec     = 0.5
sigma_pm_ra   = 0.081
sigma_pm_dec  = 0.053
sigma_acc_ra  = 0.0168
sigma_acc_dec = 0.0118
sigma_rv      = 3350.0
```

Per-star injected orbital elements live under `[mock.stars.<STAR>]`:

```toml
[mock.stars.A]
orbital_elements = [
    2092.194381543326,  # P [yr]
    0.346043688274,     # e
    2.134772053875,     # i [rad]
    2.399030601970,     # omega [rad]
    -2.382677416679,    # Omega [rad]
    -2.976490100909     # theta [rad]
]
```

The mock generator builds a Keplerian orbit from these values, evaluates it at
the reference epoch, applies the IMBH offset relative to the cluster center,
and draws noisy observations using the sigma values above.

Set `enabled = false` to fit the real `octo_utils.stars`.

### `[sampling]`

Controls the Pigeons run.

```toml
[sampling]
n_rounds             = 16
n_chains             = 192
n_chains_variational = 192
checkpoint           = true
```

`n_rounds` is the total `2^n_rounds` iterations. For Slurm runs,
`n_chains` and `julia_threads` should generally match the allocated CPU count.

### `[restart]`

Leave `job_id` empty for a fresh run.

```toml
[restart]
job_id = ""
```

To resume a checkpointed run, set `job_id` to the previous Slurm job ID and set
`[sampling].n_rounds` to the desired total round count. `submit_job.jl` searches
`mock_results/run_outputs/` for the matching `*_pt_location.txt`, verifies the
Pigeons checkpoint folder, and passes `--resume <path>` to the fitting script.

Requires the previous run to have completed with `checkpoint = true`.

### `[slurm]`

Example:

```toml
[slurm]
account       = "account"
job_name      = "mock_MAP"
nodes         = 1
cpus_per_task = 192
mem_per_cpu   = "3G"
time          = "18:00:00"
julia_module  = "julia/1.11.3"
julia_threads = 192
mail_type     = "ALL"
mail_user     = "you@institution.ca"
```

### `[paths]`

Example:

```toml
[paths]
project    = "/home/kenzhayd/projects/def-vhenault/kenzhayd/octoIMBH_env"
output_dir = "../mock_results/run_outputs"
log_dir    = "../mock_results/logs"
```

---

## Submitting a Mock Fit

At minimum, check `[stars]`, `[mock]`, `[sampling]`, `[slurm]`, and `[paths]`.

2. Preview the generated Slurm script.

```bash
cd launch_scripts/
julia submit_job.jl ../configs/my_mock_run.toml --dry-run
```

3. Submit.

```bash
julia submit_job.jl ../configs/my_mock_run.toml
```

4. Monitor.

```bash
squeue -u $USER
tail -f ../mock_results/logs/mock_MAP_<JOBID>.out
tail -f ../mock_results/logs/output_<JOBID>.log
```

The fitting script saves the chain, writes a run summary, records the Pigeons
checkpoint location, and then calls `plot_chain.jl`.

---

## Outputs

Outputs land in `[paths].output_dir`, usually `mock_results/run_outputs/`.

| File | Contents |
|------|----------|
| `*_summary.md` | Run metadata, priors, full TOML config, and checkpoint path |
| `*_chain.fits` | Full posterior chain, loadable with `Octofitter.loadchain` |
| `*_pt_location.txt` | Pigeons checkpoint folder used for restarts |
| `*_corner.png` | Corner plot of system-level parameters |
| `*_orbit_panels.png` | Sky-plane orbit panels per star and combined |
| `*_posteriors.png` | Marginal posterior histograms |
| `*_posterior_stats.txt` | Posterior medians, credible intervals, and diagnostics |
| `*_plausibility.png` | Pericenter distance/speed checks with reference radii |
| `*_phase_accel.png` | True anomaly and acceleration-alignment diagnostics |
| `*_rv_check.png` | RV posterior prediction vs measurement for RV stars |
| `*_accel_check.png` | Acceleration posterior predictive check when acceleration is used |
| `*_imbh_position.png` | IMBH position posterior density |
| `*_orbits_3d.mp4` | 3D orbit animation |





