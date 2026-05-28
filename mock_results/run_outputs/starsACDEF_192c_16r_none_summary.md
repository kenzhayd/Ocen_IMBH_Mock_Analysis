# Run Summary

- **Date:** 2026-05-28T15:58:30.720
- **Slurm Job ID:** none
- **Stars:** A, C, D, E, F
- **Reference epoch:** 55197.0 MJD (2010.0 yr)
- **Config file:** C:\Users\macke\Clusters\Ocen_IMBH_Mock_Analysis\configs\mock_64000.toml
- **Run type:** fresh

## Sampling Parameters

| Parameter | Value |
|---|---|
| n_rounds | 16 |
| n_chains | 192 |
| n_chains_variational | 192 |
| checkpoint | true |

## System Priors

| Parameter | Prior |
|---|---|
| offsetx | Uniform(-3000, 3000) |
| M | Uniform(100, 100000) |
| offsety | Uniform(-3000, 3000) |
| plx | truncated(Normal(0.19, 0.004), lower=0) |
| z_prior | Normal(0, 4558.0) AU |

## Companion Priors (defaults)

| Parameter | Prior |
|---|---|
| theta | UniformCircular() |
| e | Uniform(0.0, 0.99) |
| P | Uniform(10, 2_000_000) |
| omega | UniformCircular() |
| Omega | UniformCircular() |
| i | Sine() |

## Full Configuration

```toml
# Mock Data Configuration (See end)

# Run locally: julia --project=C:\Users\macke\Clusters\Octofitter_imbh.jl launch_scripts/octo_orbit_direct_likelihoods.jl configs/mock_default.toml

# --- Run metadata (written into the run summary) ---
[meta]
system_name = "Mock_Omega_Cen"
description = "Mock dataset 5-star direct likelihood fit"

# --- Star selection ---
# Available: A, C, D, E, F,  
[stars]
selected = ["A", "C", "D", "E", "F"]

# --- Reference epoch ---
# Reference epoch at which position offsets, PM, and acceleration are evaluated.
# Specified in decimal years; converted internally to MJD via Octofitter.years2mjd().
# 55197.0 # 2010  from runID 10836842
[epoch]
year = 2010.0


# === Priors =============================================================
# Prior strings are parsed at runtime into Distributions.jl objects.
# Supported forms:
#   "Uniform(lo, hi)"
#   "Normal(mu, sigma)"
#   "truncated(Normal(mu, sigma), lower=L)"
#   "truncated(Normal(mu, sigma), lower=L, upper=U)"
#   "Sine()"
#   "UniformCircular()"
# Underscores in numbers (e.g. 2_000_000) are allowed inside strings.
# ========================================================================

# System-level priors (shared across all companions)
[priors.system]
plx     = "truncated(Normal(0.19, 0.004), lower=0)"   # Parallax [mas]
M       = "Uniform(100, 100000)"                        # IMBH mass [solar masses]
offsetx = "Uniform(-3000, 3000)"                        # IMBH RA offset from assumed center [mas]; ±3" covers Haberle+2024 MCMC centre (0.77" NE of AvdM10)
offsety = "Uniform(-3000, 3000)"                        # IMBH Dec offset from assumed center [mas]

# Default companion (per-star) priors — applied to every star unless overridden below
[priors.companion_defaults]
P     = "Uniform(10, 2_000_000)"    # Orbital period [yr]
e     = "Uniform(0.0, 0.99)"        # Eccentricity
i     = "Sine()"                     # Inclination [rad]
omega = "UniformCircular()"          # Argument of periastron [rad]
Omega = "UniformCircular()"          # Longitude of ascending node [rad]
theta = "UniformCircular()"          # Mean anomaly at reference epoch [rad]

# --- Data selection ---
# Which observation types to include per star.  Defaults apply to all stars.
# Set to false to exclude a data type for a specific star.
# "radial_velocity" is only used when the star has RV data in octo_utils.jl.
[data.defaults]
position        = true
proper_motion   = true
acceleration    = false
radial_velocity = true
escape_velocity = false     # Häberle-style piecewise escape velocity constraint

[data.overrides.A]
radial_velocity = false

[data.overrides.C]
radial_velocity = false

[data.overrides.D]
radial_velocity = false

# Line-of-sight (z) prior — constrains each star's LOS offset from the IMBH.
# sigma_z_au is the width of a Normal(0, σ) prior in AU.
# 0.0221 pc = 4558 AU (one-dimensional positional standard deviation from Haberle)
# Omega Cen core radius ≈ 4.1 pc ≈ 845,000 AU; half-light radius ≈ 7.9 pc ≈ 1,629,000 AU.
[data.z_prior]
sigma_z_au = 4558.0            

[restart]
job_id = ""

# Per-star overrides (uncomment to disable specific data for a star):
# [data.overrides.A]
# acceleration = false

# --- Sampling ---
[sampling]
n_rounds             = 16
n_chains             = 192
n_chains_variational = 192
checkpoint           = true


# --- Slurm / HPC ---
[slurm]
account       = "def-vhenault"
job_name      = "mock_test_1"
nodes         = 1
cpus_per_task = 192
mem_per_cpu   = "3G"
time          = "16:00:00"
julia_module  = "julia/1.11.3"
julia_threads = 192
mail_type     = "ALL"
mail_user     = "Mackenzie.hayduk@smu.ca"

# --- Paths ---
# Relative paths are resolved from the directory containing this config file.
[paths]
project    = "/home/kenzhayd/projects/def-vhenault/kenzhayd/octoIMBH_env"    # --project= argument for Julia
output_dir = "../mock_results/run_outputs"        # chain files, plots, summaries
log_dir    = "../mock_results/logs"               # Slurm stdout/stderr and tee logs


# ========== MOCK DATA CONFIGURATION ==========

[mock]
enabled = true

# Mock central IMBH parameters 
M_IMBH       = 64000.021504    # Solar masses
plx     = 0.189368        # parallax [mas]
offsetx = -134.200734
offsety = 693.058436

# Observation uncertainties (σ values) for mock data
sigma_ra = 0.5       # Position RA uncertainty [mas]
sigma_dec = 0.5      # Position Dec uncertainty [mas]
# Using average RA uncertainties (Haberle): 0.038, 0.182, 0.127, 0.082, 0.025, 0.017, 0.098
sigma_pm_ra = 0.081      # Proper motion RA uncertainty [mas/yr]
# Using average Dec uncertainties (Haberle): 0.055, 0.081, 0.056, 0.061, 0.037, 0.016, 0.062
sigma_pm_dec = 0.053     # Proper motion Dec uncertainty [mas/yr]
# Using average RA Acc uncertainties (Haberle): 0.0083, 0.0239, 0.0333, 0.0177, 0.0042, 0.0038, 0.0267
sigma_acc_ra = 0.0168   # Acceleration RA uncertainty [mas/yr²]
# Using average Dec uncertainties (Haberle): 0.0098, 0.0157, 0.0123, 0.0162, 0.0075, 0.0038, 0.0170
sigma_acc_dec = 0.0118   # Acceleration Dec uncertainty [mas/yr²]
# Average of measured uncertainties 4000 m/s (F) and 2700 m/s (E)
sigma_rv = 3350.0         # Radial velocity uncertainty [m/s]


# Known orbital parameters for all mock stars
# Parameter sets are based on median posterior values from real Ocen_IMBH_analysis fits
# Ocen_IMBH_analysis: results/run_outputs/starsACDEF_192c_18r_cont_10836842_posterior_stat.txt 

# Fitting should recover these values 
# Sample index: 238052
# Target central mass: M ≈ 64000

[mock.stars]

[mock.stars.A]
orbital_elements = [6507.401, 0.41172, 2.17287, 2.56203, -2.27967, -607968.801]

[mock.stars.C]
orbital_elements = [3282.752, 0.37358, 1.15308, -1.92715, -2.87204, -105513.089]

[mock.stars.D]
orbital_elements = [8909.840, 0.08605, 1.43622, -2.58309, -1.03684, -779130.997]

[mock.stars.E]
orbital_elements = [9024.824, 0.86753, 1.23807, 2.84113, 2.32978, -1017978.028]

[mock.stars.F]
orbital_elements = [11784.581, 0.06472, 0.23269, -1.68262, -0.24758, -927902.928]




```
