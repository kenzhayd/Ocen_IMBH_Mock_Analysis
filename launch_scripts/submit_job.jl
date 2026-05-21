#!/usr/bin/env julia
"""
    submit_job.jl — Generate and submit a Slurm job from a TOML config

Usage:
    julia submit_job.jl configs/my_run.toml              # generate + submit
    julia submit_job.jl configs/my_run.toml --dry-run     # generate only (inspect before submitting)

To resume a previous run for N more rounds:
  1. Set [restart].job_id to the Slurm job ID of the run to resume.
  2. Set [sampling].n_rounds to the desired TOTAL number of rounds.
  3. Ensure checkpoint = true so the resumed run also checkpoints.
  The PT exec folder is looked up from the *_<job_id>_pt_location.txt file
  in the output directory and baked into the generated Slurm script.

The generated Slurm script is saved to the log directory for reproducibility.
"""

using TOML
using Dates

# ── Parse arguments ─────────────────────────────────────────────────────

if isempty(ARGS)
    println(stderr, "Usage: julia submit_job.jl <config.toml> [--dry-run]")
    exit(1)
end

config_path = ARGS[1]
dry_run = "--dry-run" in ARGS

isfile(config_path) || error("Config file not found: $config_path")
cfg = TOML.parsefile(config_path)

# ── Extract sections ────────────────────────────────────────────────────

slurm   = cfg["slurm"]
paths   = cfg["paths"]
restart = get(cfg, "restart", Dict())

# Resolve paths relative to the config file's directory
config_dir = dirname(abspath(config_path))
abs_config = abspath(config_path)

log_dir     = isabspath(paths["log_dir"])     ? paths["log_dir"]     : joinpath(config_dir, paths["log_dir"])
output_dir  = isabspath(paths["output_dir"])  ? paths["output_dir"]  : joinpath(config_dir, paths["output_dir"])
project_dir = isabspath(paths["project"])     ? paths["project"]     : joinpath(config_dir, paths["project"])

# ── Optional resume: resolve PT exec folder from job_id ─────────────────
# The fitting script writes a *_<job_id>_pt_location.txt file in output_dir
# at the end of each checkpointed run.  Supply that job_id here to resume.
resume_arg = ""
job_id_str = get(restart, "job_id", "")
if job_id_str isa String && !isempty(job_id_str)
    mkpath(output_dir)
    candidates = filter(
        f -> endswith(f, "_pt_location.txt") && occursin("_$(job_id_str)_", f),
        readdir(output_dir; join=true)
    )
    isempty(candidates) && error(
        "No pt_location file found for job_id=$(job_id_str) in $(output_dir). " *
        "Ensure the previous run completed with checkpoint=true."
    )
    length(candidates) > 1 && @warn(
        "Multiple pt_location files match job_id=$(job_id_str); using: $(candidates[1])"
    )
    pt_exec_folder = strip(read(candidates[1], String))
    isdir(pt_exec_folder) || error("PT exec folder does not exist: $pt_exec_folder")
    resume_arg = " --resume $(pt_exec_folder)"
    println("Resume job_id=$(job_id_str) → PT folder: $pt_exec_folder")
end

# The fitting script lives next to this launcher
fitting_script = joinpath(@__DIR__, "octo_orbit_direct_likelihoods.jl")

# ── Generate Slurm script ──────────────────────────────────────────────

job_name = slurm["job_name"]

script = """
#!/bin/bash
#SBATCH --account=$(slurm["account"])
#SBATCH --job-name=$(job_name)
#SBATCH --nodes=$(slurm["nodes"])
#SBATCH --cpus-per-task=$(slurm["cpus_per_task"])
#SBATCH --mem-per-cpu=$(slurm["mem_per_cpu"])
#SBATCH --time=$(slurm["time"])
#SBATCH --output=$(log_dir)/$(job_name)_%j.out
#SBATCH --error=$(log_dir)/$(job_name)_%j.err
#SBATCH --mail-type=$(slurm["mail_type"])
#SBATCH --mail-user=$(slurm["mail_user"])

mkdir -p $(log_dir)
mkdir -p $(output_dir)

module load $(slurm["julia_module"])

export JULIA_CONDAPKG_BACKEND=Null

julia --project=$(project_dir) -e 'using Pkg; Pkg.instantiate(); Pkg.add(["CairoMakie", "PairPlots", "Distributions", "Unitful", "UnitfulAstro", "Pigeons"])' 2>&1 | tee $(log_dir)/instantiate_\${SLURM_JOB_ID}.log

julia --project=$(project_dir) -t $(slurm["julia_threads"]) \\
    $(fitting_script) \\
    $(abs_config)$(resume_arg) 2>&1 | tee $(log_dir)/output_\${SLURM_JOB_ID}.log
"""

# ── Write the script ───────────────────────────────────────────────────

mkpath(log_dir)
timestamp = Dates.format(now(), "yyyymmdd_HHMMSS")
script_path = joinpath(log_dir, "job_$(timestamp).sh")
write(script_path, script)
println("Generated Slurm script: $script_path")

# ── Submit (or not) ────────────────────────────────────────────────────

if dry_run
    println("Dry run — not submitting. Inspect the script above, then run:")
    println("  sbatch $script_path")
else
    run(`sbatch $script_path`)
end
