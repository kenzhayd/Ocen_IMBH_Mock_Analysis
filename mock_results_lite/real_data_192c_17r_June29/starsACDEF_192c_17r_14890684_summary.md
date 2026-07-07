# Run Summary

- **Date:** 2026-06-29T02:37:43.223
- **Slurm Job ID:** 14890684
- **Stars:** A, C, D, E, F
- **Reference epoch:** 55197.0 MJD (2010.0 yr)
- **Config file:** /home/kenzhayd/projects/def-vhenault/kenzhayd/Ocen_IMBH_Mock_Analysis/configs/10836842_rerun.toml
- **Run type:** fresh

## Sampling Parameters

| Parameter | Value |
|---|---|
| n_rounds | 17 |
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
# ============================================================
# Omega Centauri IMBH Orbit Fitting — Run Configuration 
# Identical to run ID: 10836842 - Primary Fit
# ============================================================
#
[meta]
system_name = "Omega_Cen"
description = "Rerun of primary 5-star direct likelihood fit (no accel, z-prior)"

[stars]
selected = ["A", "C", "D", "E", "F"]

[epoch]
year = 2010.0

[priors.system]
plx     = "truncated(Normal(0.19, 0.004), lower=0)"
M       = "Uniform(100, 100000)"
offsetx = "Uniform(-3000, 3000)"
offsety = "Uniform(-3000, 3000)"

[priors.companion_defaults]
P     = "Uniform(10, 2_000_000)"
e     = "Uniform(0.0, 0.99)"
i     = "Sine()"
omega = "UniformCircular()"
Omega = "UniformCircular()"
theta = "UniformCircular()"

[data.defaults]
position        = true
proper_motion   = true
acceleration    = false
radial_velocity = true

z_prior         = true
escape_velocity = false

[data.z_prior]
sigma_z_au = 4558

[restart]
job_id = ""

[sampling]
n_rounds             = 17
n_chains             = 192
n_chains_variational = 192
checkpoint           = true

[slurm]
account       = "def-vhenault"
job_name      = "octo_imbh"
nodes         = 1
cpus_per_task = 192
mem_per_cpu   = "3G"
time          = "23:59:00"
julia_module  = "julia/1.11.3"
julia_threads = 192
mail_type     = "ALL"
mail_user     = "mackenzie.hayduk@smu.ca"

[paths]
project    = "/home/kenzhayd/projects/def-vhenault/kenzhayd/octoIMBH_env"    # --project= argument for Julia
output_dir = "../results/run_outputs"        # chain files, plots, summaries
log_dir    = "../results/logs"               # Slurm stdout/stderr and tee logs


```

## Sampling Result


## Pigeons PT Checkpoint

The PT exec folder contains intermediate checkpoint files.
It can be deleted once the chain FITS file has been verified.

| | Path |
|---|---|
| PT exec folder | `/home/kenzhayd/projects/def-vhenault/kenzhayd/Ocen_IMBH_Mock_Analysis/launch_scripts/results/all/2026-06-29-02-38-26-h9tJL4Iu` |
| PT location file | `/home/kenzhayd/projects/def-vhenault/kenzhayd/Ocen_IMBH_Mock_Analysis/configs/../results/run_outputs/starsACDEF_192c_17r_14890684_pt_location.txt` |
