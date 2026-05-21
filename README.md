# Ocen IMBH Analysis

Scripts for fitting orbits of high-velocity stars around an intermediate-mass black hole (IMBH) candidate in ω Centauri. Orbit fitting uses [`Octofitter_imbh.jl`](../Octofitter_imbh.jl), a development fork of [Octofitter](https://github.com/sefffal/Octofitter.jl), sampling with Pigeons (parallel tempering HMC/NUTS).

---

## Repository Structure

```
Ocen_IMBH_analysis/
├── configs/
│   └── default.toml                          # Reference run configuration
├── launch_scripts/
│   ├── octo_utils.jl                         # Star data, RV constants, observation builder
│   ├── octo_orbit_direct_likelihoods.jl      # Main fitting script (v8 API)
│   ├── parse_config.jl                       # TOML loader, prior parser, data flag helpers
│   ├── plot_chain.jl                         # Post-fit plotting (posteriors, orbits, 3D animation)
│   ├── submit_job.jl                         # Slurm job generator/submitter
│   └── old_scripts/
│       ├── octo_orbit_192c_18r.jl            # Legacy fitting script (v7 API)
│       └── run_octo_192_all_MCMC_centre.sh   # Legacy Slurm launcher
├── results/
│   ├── run_outputs/                          # Chain files, plots, summaries
│   │   └── run_plot_chain.sh                 # Slurm job to regenerate plots from a chain
│   └── logs/                                 # Slurm stdout/stderr and generated job scripts
└── tests/
    ├── run_tests.sh                          # Slurm batch script to run all tests
    ├── test_likelihoods.jl                   # Unit tests for PM/accel likelihoods
    ├── test_rv_likelihood.jl                 # Unit tests for RV likelihood (PlanetRelativeRVObs)
    ├── test_ad_compatibility.jl              # ForwardDiff gradient verification
    └── test_small_fit.jl                     # End-to-end 2-star sampling test
```

---

## Prerequisites

- Julia ≥ 1.10 with the `Octofitter_imbh.jl` package (the sibling directory `../Octofitter_imbh.jl`)
- On DAC HPC clusters: `module load julia/1.10.10`

All scripts are run with `--project=path/to/Octofitter_imbh.jl` so that the local fork is used rather than any registered Octofitter version.

The fitting script auto-installs `OctofitterRadialVelocity` (a sub-package inside `Octofitter_imbh.jl`) if it is not already present in the project environment.

---

## Configuration File

All run parameters live in a single TOML file. The reference configuration is [`configs/default.toml`](configs/default.toml). Copy and edit it to create a new run variant:

```bash
cp configs/default.toml configs/my_run.toml
# edit configs/my_run.toml
```

### Key sections

#### `[meta]`
Human-readable labels written into the run summary.

```toml
[meta]
system_name = "Omega_Cen"
description = "5-star direct likelihood fit"
```

#### `[stars]`
Which stars to include. Available labels: `A`, `B`, `C`, `D`, `E`, `F`, `G`. Stars B and G have lower-quality data. Stars E and F have radial velocity measurements.

```toml
[stars]
selected = ["A", "C", "D", "E", "F"]
```

#### `[epoch]`
Reference epoch at which position, proper motion, and acceleration are evaluated (decimal years; converted to MJD internally).

```toml
[epoch]
year = 2010.0
```

#### `[priors.system]`
System-level parameters shared across all stars.

| Parameter | Meaning | Units |
|-----------|---------|-------|
| `plx` | Parallax | mas |
| `M` | IMBH mass | M☉ |
| `offsetx` | IMBH RA offset from assumed centre | mas |
| `offsety` | IMBH Dec offset from assumed centre | mas |

```toml
[priors.system]
plx     = "truncated(Normal(0.19, 0.004), lower=0)"
M       = "Uniform(100, 100000)"
offsetx = "Uniform(-3000, 3000)"
offsety = "Uniform(-3000, 3000)"
```

#### `[priors.companion_defaults]`
Default orbital priors applied to every star.

| Parameter | Meaning |
|-----------|---------|
| `P` | Period [yr] |
| `e` | Eccentricity |
| `i` | Inclination (use `"Sine()"`) |
| `omega` | Argument of periastron (use `"UniformCircular()"`) |
| `Omega` | Longitude of ascending node (use `"UniformCircular()"`) |
| `theta` | Mean anomaly at epoch (use `"UniformCircular()"`) |

#### `[priors.overrides.<STAR>]`  *(optional)*
Override defaults for individual stars. Only list parameters that differ.

```toml
[priors.overrides.B]
P = "Uniform(1, 2_000_000)"
```

#### Prior specification syntax

Prior strings are parsed at runtime. Supported forms:

```
"Uniform(lo, hi)"
"Normal(mu, sigma)"
"truncated(Normal(mu, sigma), lower=L)"
"truncated(Normal(mu, sigma), lower=L, upper=U)"
"Sine()"
"UniformCircular()"
```

Underscores in numbers are allowed (`"Uniform(10, 800_000)"`).

#### `[data.defaults]` and `[data.overrides.<STAR>]`
Controls which observation types enter the likelihood per star. Defaults apply to all stars; per-star overrides can disable specific data types. Radial velocity is only used for stars that have RV data in `octo_utils.jl` (currently E and F).

```toml
[data.defaults]
position        = true
proper_motion   = true
acceleration    = true
radial_velocity = true
z_prior         = true      # LOS position prior (cluster membership)
escape_velocity = false     # Häberle-style piecewise escape velocity constraint

[data.z_prior]
sigma_z_au = 845_000        # Normal(0, σ) prior on LOS offset [AU]

# Per-star overrides (uncomment to disable specific data for a star):
# [data.overrides.A]
# acceleration = false
```

The `z_prior` adds a `Normal(0, sigma_z_au)` prior on each star's line-of-sight separation from the IMBH (in AU), constraining stars to lie within the cluster extent. The default `sigma_z_au = 845,000` corresponds to the core radius of ω Cen (~4.1 pc).

The `escape_velocity` flag adds a `PlanetEscVelObs` likelihood following Häberle et al. (2024). It applies no penalty when the modelled total escape velocity at the star's 3D position exceeds the observed 2D speed (the star is bound), and a Gaussian penalty in the velocity excess otherwise:

```
ln ℒ = 0                                   if v₂D ≤ v_esc,total
ln ℒ = −½ ((v₂D − v_esc,total) / σ_v)²   if v₂D > v_esc,total
```

where `v_esc,total = √(v_esc,IMBH² + v_esc,cluster²)`, `v_esc,IMBH = √(2GM/r₃D)`, and `r₃D` is the 3D separation computed from the orbital solution. The observed 2D speed `v₂D` and its uncertainty `σ_v` come from `star.v2D` / `star.v2D_err` in `octo_utils.jl`. The cluster escape velocity is fixed at `v_esc_cluster = 62.0 km/s` (Häberle et al. 2024). This likelihood is disabled by default; enable it with `escape_velocity = true`.

#### `[sampling]`
Controls Pigeons parallel-tempering run.

```toml
[sampling]
n_rounds             = 15     # 2^n_rounds total iterations
n_chains             = 192    # parallel chains (should match CPU count)
n_chains_variational = 192
checkpoint           = true   # must be true to enable restarts
```

#### `[restart]` *(optional)*
Resumes a previous run for additional rounds. Leave `job_id` empty (or omit the section entirely) for a fresh run.

```toml
[restart]
job_id = ""   # Slurm job ID of the run to resume (e.g. "9781546")
```

When `job_id` is set:
- `submit_job.jl` looks up `results/run_outputs/*_<job_id>_pt_location.txt` to find the Pigeons checkpoint directory and bakes `--resume <path>` into the generated Slurm script.
- `n_rounds` in `[sampling]` is the **total** round target (e.g. set to `20` to add 5 rounds to a previous 15-round run). The script errors if `n_rounds` ≤ rounds already completed.
- The resumed run also writes a new `*_pt_location.txt` and appends its checkpoint path to the run summary, so the chain can be extended again if needed.

#### `[slurm]`
Slurm resource requests. Adjust for the target cluster.

```toml
[slurm]
account       = "def-vhenault"
job_name      = "octo_imbh"
nodes         = 1
cpus_per_task = 192
mem_per_cpu   = "3G"
time          = "15:00:00"
julia_module  = "julia/1.10.10"
julia_threads = 192
mail_type     = "ALL"
mail_user     = "you@institution.ca"
```

#### `[paths]`
Paths are resolved relative to the config file's directory.

```toml
[paths]
project    = "../../Octofitter_imbh.jl"     # --project= for Julia
output_dir = "../results/run_outputs"        # chain files, plots, summaries
log_dir    = "../results/logs"               # Slurm stdout/stderr and tee logs
```

---

## Running Tests

The test scripts can be run individually or via the Slurm batch script. Since `test_small_fit.jl` runs a short MCMC it benefits from multiple CPUs and is best submitted as a cluster job.

### Run on HPC (recommended)

```bash
cd tests/
sbatch run_tests.sh
```

Before submitting, edit `run_tests.sh` to set the correct `OCTOFITTER_PROJECT` path:

```bash
OCTOFITTER_PROJECT=/path/to/Octofitter_imbh.jl
```

Output is written to `tests/logs/octo_tests_<JOBID>.out` and `.err`.

### Run locally (fast tests only)

`test_likelihoods.jl`, `test_ad_compatibility.jl`, and `test_rv_likelihood.jl` are lightweight and finish in under a minute:

```bash
julia --project=../Octofitter_imbh.jl tests/test_likelihoods.jl
julia --project=../Octofitter_imbh.jl tests/test_ad_compatibility.jl
julia --project=../Octofitter_imbh.jl tests/test_rv_likelihood.jl
```

`test_small_fit.jl` runs a real 2-star fit and requires multiple threads:

```bash
julia --project=../Octofitter_imbh.jl -t 4 tests/test_small_fit.jl
```

### What the tests check

| Script | What it tests |
|--------|--------------|
| `test_likelihoods.jl` | `PlanetPMObs` and `PlanetAccelObs` evaluate correctly on a known orbit; position/PM/acceleration recover injected values |
| `test_rv_likelihood.jl` | `PlanetRelativeRVObs` construction, combined pos+PM+accel+RV model compilation, ForwardDiff gradient through RV likelihood |
| `test_ad_compatibility.jl` | Model compiles as type-stable; ForwardDiff can compute gradients without errors |
| `test_small_fit.jl` | End-to-end: 2-star MCMC sampling completes; posterior plots are generated |

---

## Submitting a Production Fit

### Fresh run

#### 1. Create (or choose) a config file

```bash
cp configs/default.toml configs/my_run.toml
# edit configs/my_run.toml as needed
```

Adjust at minimum: `[stars] selected`, `[slurm] account`, `[slurm] mail_user`, and `[paths] project`. Set `checkpoint = true` in `[sampling]` to enable future restarts.

#### 2. Preview the generated Slurm script (dry run)

```bash
cd launch_scripts/
julia submit_job.jl ../configs/mock_default.toml --dry-run
```

This writes a timestamped `.sh` file to the log directory and prints the `sbatch` command, but does **not** submit.

#### 3. Submit

```bash
julia submit_job.jl ../configs/my_run.toml
```

`submit_job.jl` parses the config, generates a Slurm batch script in `<log_dir>/job_<timestamp>.sh`, and calls `sbatch`. The job runs `octo_orbit_direct_likelihoods.jl` with the config as its argument; `plot_chain.jl` is called automatically at the end to generate all plots and summaries.

#### 4. Monitor

```bash
squeue -u $USER
tail -f results/logs/octo_imbh_<JOBID>.out
```

---

### Restarting a run for additional rounds

Requires the original run to have completed with `checkpoint = true`.

#### 1. Find the job ID

Look for `*_pt_location.txt` in `results/run_outputs/` or check the original run summary — the job ID is the numeric suffix in every output filename (e.g. `starsACDEF_192c_15r_9781546_chain.fits` → job ID `9781546`).

#### 2. Edit the config

```toml
[restart]
job_id = "9781546"   # job ID of the run to resume

[sampling]
n_rounds  = 20       # TOTAL target (e.g. 20 adds 5 rounds to a 15-round run)
checkpoint = true
```

#### 3. Submit

```bash
julia submit_job.jl ../configs/my_run.toml
```

`submit_job.jl` reads the `pt_location.txt` file, verifies the checkpoint directory exists, and bakes `--resume <path>` into the generated Slurm script. The resumed run produces outputs with a `_cont` suffix (e.g. `starsACDEF_192c_20r_cont_<NEWJOBID>_chain.fits`) and writes its own `*_pt_location.txt` so the chain can be extended again.

---

### 5. Outputs

All outputs land in `<output_dir>` (default: `results/run_outputs/` relative to the config file):

| File | Contents |
|------|---------|
| `*_summary.md` | Run metadata, priors, full config, and PT checkpoint path |
| `*_chain.fits` | Full posterior chain (load with `Octofitter.loadchain`) |
| `*_pt_location.txt` | Path to the Pigeons PT checkpoint directory (used for restarts; can be deleted once the chain FITS is verified) |
| `*_corner.png` | Corner plot of system-level parameters |
| `*_orbit_panels.png` | Sky-plane orbit panels per star + combined |
| `*_posteriors.png` | Marginal posterior histograms |
| `*_posterior_stats.txt` | Posterior summaries (median + 68% CI) including physical diagnostics |
| `*_plausibility.png` | Pericenter distance and speed distributions with tidal/Schwarzschild reference lines |
| `*_phase_accel.png` | True anomaly at obs epoch and acceleration-vector misalignment angle |
| `*_rv_check.png` | RV posterior prediction vs measurement (stars with RV data only) |
| `*_accel_check.png` | Acceleration posterior predictive check: predicted vs measured (accra, accdec) per star |
| `*_imbh_position.png` | IMBH position 2D density map |
| `*_orbits_3d.mp4` | 3D orbit animation (360° pan with elevation oscillation) |

---

## Regenerating Plots from an Existing Chain

To regenerate plots without re-running the fit, use the standalone Slurm script:

```bash
cd results/run_outputs/
sbatch run_plot_chain.sh <chain.fits>
```

This submits a lightweight 1-CPU job that runs `plot_chain.jl` on the specified chain file.

---

## Running the Fitting Script Directly

For local testing or interactive use, bypass `submit_job.jl` and call the fitting script directly:

```bash
# Fresh run
julia --project=../Octofitter_imbh.jl -t <threads> \
    launch_scripts/octo_orbit_direct_likelihoods.jl \
    configs/my_run.toml

# Resume from a checkpoint (n_rounds in config = total target)
julia --project=../Octofitter_imbh.jl -t <threads> \
    launch_scripts/octo_orbit_direct_likelihoods.jl \
    configs/my_run.toml --resume /path/to/pt_exec_folder
```

If no config path is given, the script falls back to `configs/default.toml`.

---

## Utility Module

[`launch_scripts/octo_utils.jl`](launch_scripts/octo_utils.jl) is a Julia module providing:

- `StarData` struct — position (RA/Dec), proper motion, acceleration, 2D velocity, and radial velocity with uncertainties for each star
- `stars` dictionary — labelled `"A"` through `"G"` with measured observables
- `build_star_observations(star, epoch_mjd; include_rv, z_prior_sigma, include_esc_vel)` — returns `(astrom_obs, pm_obs, acc_obs, rv_obs, zp_obs, ev_obs)`; `rv_obs`, `zp_obs`, and `ev_obs` are `nothing` when the star has no RV data, when `z_prior_sigma` is not supplied, or when `include_esc_vel=false` (the default)
- Cluster constants: centre coordinates (Anderson & van der Marel 2010), distance, systemic radial velocity (Baumgardt catalogue), central escape velocity `v_esc_cluster_kms = 62.0 km/s` (Häberle et al. 2024)
- Unit conversion and error propagation utilities

This module is included automatically by `octo_orbit_direct_likelihoods.jl` and all test scripts.

---

## Plot Chain Script

[`launch_scripts/plot_chain.jl`](launch_scripts/plot_chain.jl) loads a saved chain FITS file and generates:

1. **Corner plot** (`*_corner.png`) — pair plot of system-level parameters (M, plx, offsetx, offsety)
2. **Orbit panels** (`*_orbit_panels.png`) — sky-plane orbits per star with observed proper motion vectors, plus a combined panel
3. **Posterior histograms** (`*_posteriors.png`) — marginal distributions for each star's orbital elements
4. **Posterior statistics** (`*_posterior_stats.txt`) — median and 68% credible intervals for all orbital elements and physical diagnostics; printed to screen and saved to file
5. **Physical plausibility** (`*_plausibility.png`) — per-star histograms of pericenter distance (log-scale AU) and pericenter speed (km/s via vis-viva), with reference lines for the main-sequence tidal disruption radius, red-giant tidal disruption radius, and Schwarzschild radius
6. **Phase and alignment** (`*_phase_accel.png`) — per-star histograms of (1) the true anomaly at the observation epoch ν(t_obs), (2) the misalignment angle Δφ between the measured acceleration vector and the star→IMBH direction with a shaded band showing the measurement angular uncertainty σ_φ, and (3) the z-score of the measured acceleration component toward the IMBH (accounts for measurement SNR; |z| < 1 means the measurement is direction-agnostic)
7. **RV consistency** (`*_rv_check.png`, only for stars with RV data) — posterior-predicted peculiar radial velocity compared to the measured value ± 1σ; only produced for stars E and F
8. **Acceleration posterior predictive check** (`*_accel_check.png`) — for each star: (left) 2D scatter of posterior-predicted (accra, accdec) with the measured value overlaid as a cross with ±1σ error bars; (right) histogram of the per-sample 2D chi-squared residual χ²_i = (accra_meas − accra_pred)²/σ_ra² + (accdec_meas − accdec_pred)²/σ_dec², with a dashed line at χ² = 2.30 (the 68% contour of a 2-DOF chi-squared distribution). The fraction of draws inside the 1σ ellipse (f₆₈) is shown in the legend and recorded in `*_posterior_stats.txt`. A value near 0.68 indicates consistency; f₆₈ ≪ 0.68 indicates the model cannot reproduce the measured acceleration.
9. **IMBH position map** (`*_imbh_position.png`) — 2D density of the IMBH position posterior with absolute RA/Dec secondary axes
10. **3D orbit animation** (`*_orbits_3d.mp4`) — rotating view of all orbits in parsec, IMBH-centric frame, with star marker sizes scaled by the IMBH mass of each chain sample

The `*_posterior_stats.txt` file includes orbital elements, pericenter/apocenter distances, pericenter/apocenter speeds, orbital periods, tidal radii, Schwarzschild radius, true anomaly at epoch, acceleration misalignment angle Δφ, acceleration direction uncertainty σ_φ, acceleration z-score toward the IMBH, acceleration predictive check median χ² and f₆₈ per star, and (where applicable) the RV residual in units of σ.

This script is called automatically at the end of a fitting run, or can be run standalone via `run_plot_chain.sh`.

---

## Legacy Scripts

[`launch_scripts/old_scripts/`](launch_scripts/old_scripts/) contains:

- `octo_orbit_192c_18r.jl` — v7 Octofitter API with synthetic 3-epoch astrometry rather than direct PM/acceleration likelihoods
- `run_octo_192_all_MCMC_centre.sh` — corresponding Slurm launcher

These are kept for comparison only and are not actively maintained.
