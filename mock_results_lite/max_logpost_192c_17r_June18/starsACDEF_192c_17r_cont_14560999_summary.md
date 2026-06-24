# Run Summary

- **Date:** 2026-06-22T06:52:15.751
- **Slurm Job ID:** 14560999
- **Stars:** A, C, D, E, F
- **Reference epoch:** 55197.0 MJD (2010.0 yr)
- **Config file:** /home/kenzhayd/projects/def-vhenault/kenzhayd/Ocen_IMBH_Mock_Analysis/configs/mock_max_logpost_restart2.toml
- **Run type:** continuation (target: 17 total rounds)
- **Resumed from PT checkpoint:** /home/kenzhayd/projects/def-vhenault/kenzhayd/Ocen_IMBH_Mock_Analysis/launch_scripts/results/all/2026-06-21-18-00-42-GNCUY8Ep

## Sampling Parameters

| Parameter | Value |
|---|---|
| n_rounds | 17 (total) |
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
job_id = "14524687"

# Per-star overrides (uncomment to disable specific data for a star):
# [data.overrides.A]
# acceleration = false

# --- Sampling ---
[sampling]
n_rounds             = 17
n_chains             = 192
n_chains_variational = 192
checkpoint           = true


# --- Slurm / HPC ---
[slurm]
account       = "def-vhenault"
job_name      = "mock_max_logpost"
nodes         = 1
cpus_per_task = 192
mem_per_cpu   = "3G"
time          = "18:00:00"
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
          
# System parameters 
M_IMBH  = 63428.188701256993    # Solar masses
plx     = 0.194037631899    # parallax [mas]
offsetx = -136.391981772081
offsety = 772.678117346523

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


# Posterior orbital parameters from Ocen_IMBH_analysis
#Sample index 243695 of
# starsACDEF_192c_18r_cont_10836842_chain.fits
# period, eccentricity, inclination [rad], ω [rad], Ω [rad], theta [rad]


[mock.stars]

[mock.stars.A]
orbital_elements = [1742.492266784197, 0.560668674636, 2.346848256902, 3.009930277526, -2.196517158486, -2.992019583391]

[mock.stars.C]
orbital_elements = [2253.894969208533, 0.047279368601, 1.335279239253, 3.057377067795, -2.823694509911, -1.702768819719]

[mock.stars.D]
orbital_elements = [3086.704980821688, 0.084331429717, 1.474639751894, -1.532243035103, -1.001863515711, 0.888610279565]

[mock.stars.E]
orbital_elements = [7119.894533498357, 0.502670589295, 1.408540530364, 2.464969446310, 2.285994597414, 2.410224493155]

[mock.stars.F]
orbital_elements = [5405.142037710048, 0.068286787073, 0.425628982041, -1.137072201885, -0.088752602761, 1.426043923174]

```

## Sampling Result

- **Additional rounds run:** 1
- **Total rounds:** 17

## Pigeons PT Checkpoint

The PT exec folder contains intermediate checkpoint files.
It can be deleted once the chain FITS file has been verified.

| | Path |
|---|---|
| PT exec folder | `/home/kenzhayd/projects/def-vhenault/kenzhayd/Ocen_IMBH_Mock_Analysis/launch_scripts/results/all/2026-06-22-06-52-16-KWwcBQ49` |
| PT location file | `/home/kenzhayd/projects/def-vhenault/kenzhayd/Ocen_IMBH_Mock_Analysis/configs/../mock_results/run_outputs/starsACDEF_192c_17r_cont_14560999_pt_location.txt` |
