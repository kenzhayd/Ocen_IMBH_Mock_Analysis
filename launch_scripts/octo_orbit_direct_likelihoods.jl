"""
Orbital Models of High Velocity Stars in Omega Centauri
Using Octofitter — Direct PM & Acceleration Likelihoods

Uses direct single-epoch position, proper motion, and acceleration
likelihoods instead of synthetic multi-epoch astrometry.

Usage:
    # Fresh run:
    julia --project=<Octofitter_imbh.jl> -t <threads> octo_orbit_direct_likelihoods.jl <config.toml>

    # Resume from a Pigeons checkpoint (run n_rounds additional rounds):
    julia --project=<Octofitter_imbh.jl> -t <threads> octo_orbit_direct_likelihoods.jl <config.toml> --resume <pt_exec_folder>

    # Alternatively, set [restart].job_id in the config to the Slurm job ID of
    # the run to resume (submit_job.jl resolves the PT folder automatically).

If no config path is given, falls back to ../configs/default.toml.
"""

ENV["OCTOFITTERPY_AUTOLOAD_EXTENSIONS"] = "yes"

# Ensure OctofitterRadialVelocity is dev'd into the project environment
# (it is a sub-package inside the Octofitter repo, not auto-discovered).
import Pkg
let rv_pkg = "OctofitterRadialVelocity"
    deps = Pkg.dependencies()
    if !any(p -> p.second.name == rv_pkg, deps)
        rv_path = normpath(joinpath(@__DIR__, "..", "..", "Octofitter_imbh.jl", "OctofitterRadialVelocity"))
        @info "Adding $rv_pkg from $rv_path"
        Pkg.develop(path=rv_path)
    end
end

using Octofitter
using Octofitter: @variables, System
using Distributions
using Unitful
using UnitfulAstro
using LinearAlgebra
using Statistics
using Dates
using Pigeons
using OctofitterRadialVelocity

# Add the directory to LOAD_PATH
push!(LOAD_PATH, @__DIR__)
using octo_utils  # local module

# Load configuration helpers and config file
include(joinpath(@__DIR__, "parse_config.jl"))
config_path = length(ARGS) >= 1 ? ARGS[1] : joinpath(@__DIR__, "..", "configs", "default.toml")
cfg = load_config(config_path)
println("Loaded config: $config_path")

# ── Restart / resume detection ───────────────────────────────────────────
# --resume <path> CLI flag takes precedence over [restart].job_id in config.
# In both cases the PT exec folder is located from a *_pt_location.txt file
# written at the end of each run (or supplied directly via --resume <path>).
resume_pt_folder = nothing
let i = findfirst(==("--resume"), ARGS)
    if i !== nothing && i < length(ARGS)
        resume_pt_folder = ARGS[i+1]
    end
end
if resume_pt_folder === nothing
    restart_cfg = get(cfg, "restart", Dict())
    job_id_str  = get(restart_cfg, "job_id", "")
    if job_id_str isa String && !isempty(job_id_str)
        # Resolve output_dir early so we can search for the pt_location file.
        paths_cfg_early = cfg["paths"]
        output_dir_early = isabspath(paths_cfg_early["output_dir"]) ?
            paths_cfg_early["output_dir"] :
            joinpath(dirname(abspath(config_path)), paths_cfg_early["output_dir"])
        candidates = filter(
            f -> endswith(f, "_pt_location.txt") && occursin("_$(job_id_str)_", f),
            readdir(output_dir_early; join=true)
        )
        if isempty(candidates)
            error("No pt_location file found for job_id=$(job_id_str) in $(output_dir_early). " *
                  "Ensure the previous run completed with checkpoint=true and wrote a *_pt_location.txt file.")
        elseif length(candidates) > 1
            @warn "Multiple pt_location files found for job_id=$(job_id_str); using the first: $(candidates[1])"
        end
        resume_pt_folder = strip(read(candidates[1], String))
    end
end
is_resume = resume_pt_folder !== nothing
if is_resume
    println("Resume mode: will load PT checkpoint from: $resume_pt_folder")
    isdir(resume_pt_folder) || error("PT checkpoint folder not found: $resume_pt_folder")
end

# === 1. Select stars and time config ===
star_names = cfg["stars"]["selected"]
epoch_mjd  = get_epoch_mjd(cfg)
epoch_year = cfg["epoch"]["year"]

# ========== START MOCK MODIFICATIONS ==========
# Determine whether to use mock data or real data
is_mock_enabled = haskey(cfg, "mock") && get(cfg["mock"], "enabled", false)
if is_mock_enabled
    println("Mock data fitting: generating synthetic data from config parameterset")
end


# === 2. Build observation objects for each star ===
astrom_obs = Dict{String, Any}()
pm_obs     = Dict{String, Any}()
acc_obs    = Dict{String, Any}()
rv_obs     = Dict{String, Any}()
zp_obs     = Dict{String, Any}()
ev_obs     = Dict{String, Any}()

z_prior_sigma = get_z_prior_sigma(cfg)

for name in star_names
    # ========== LOAD DATA ==========
    # Either from mock config (if enabled) or from real data dictionary
    if is_mock_enabled
        # Generate mock StarData from orbital parameters in config
        mock_cfg = cfg["mock"]
        star_cfg = mock_cfg["stars"][name]
        a, e, i, ω, Ω, tp = star_cfg["orbital_elements"]
        
        star = octo_utils.stardata_struct(
            name;
            a=a, e=e, i=i, ω=ω, Ω=Ω, tp=tp,
            M=mock_cfg["M_IMBH"],
            plx=mock_cfg["plx"],
            t_ref=epoch_mjd,
            epoch=epoch_mjd,
            sigma_ra=mock_cfg["sigma_ra"],
            sigma_dec=mock_cfg["sigma_dec"],
            sigma_pm_ra=mock_cfg["sigma_pm_ra"],
            sigma_pm_dec=mock_cfg["sigma_pm_dec"],
            sigma_acc_ra=mock_cfg["sigma_acc_ra"],
            sigma_acc_dec=mock_cfg["sigma_acc_dec"],
            sigma_rv=mock_cfg["sigma_rv"]
        )
    else
        # Load real star data from dictionary
        star = octo_utils.stars[name]
    end
    
    
    include_rv  = get_data_flag(cfg, name, "radial_velocity")
    include_ev  = get_data_flag(cfg, name, "escape_velocity")
    include_acc = get_data_flag(cfg, name, "acceleration")
    
    
    
    # use mock data (with noise) or real data (direct values)
    if is_mock_enabled
        a, p, ac, r, zp, ev = octo_utils.build_mock_observations(star, epoch_mjd;
                            include_rv, z_prior_sigma, include_esc_vel=include_ev, include_acc=include_acc)
    else
        a, p, ac, r, zp, ev = octo_utils.build_star_observations(star, epoch_mjd;
                            include_rv, z_prior_sigma, include_esc_vel=include_ev,
                            include_acc=include_acc)
    end 


    astrom_obs[name] = a
    pm_obs[name]     = p
    acc_obs[name]    = ac
    rv_obs[name]     = r
    zp_obs[name]     = zp
    ev_obs[name]     = ev
end

# ========== END MOCK MODIFICATIONS ========== 


# === 3. Define companions ===
companions = Planet[]
for name in star_names
    # Parse priors from config (with per-star overrides)
    P_prior = parse_prior(get_companion_prior(cfg, name, "P"))
    e_prior = parse_prior(get_companion_prior(cfg, name, "e"))
    i_prior = parse_prior(get_companion_prior(cfg, name, "i"))
    ω_prior = parse_prior(get_companion_prior(cfg, name, "omega"))
    Ω_prior = parse_prior(get_companion_prior(cfg, name, "Omega"))
    θ_prior = parse_prior(get_companion_prior(cfg, name, "theta"))

    # Build observation list based on config data flags
    obs_list = Any[]
    if get_data_flag(cfg, name, "position")
        use_oneil = get(get(get(cfg, "data", Dict()), "defaults", Dict()), "position_oneil", true)
        if use_oneil
            push!(obs_list, ObsPriorAstromONeil2019(astrom_obs[name]))
        else
            push!(obs_list, astrom_obs[name])
        end
    end
    if get_data_flag(cfg, name, "proper_motion")
        push!(obs_list, pm_obs[name])
    end
    if acc_obs[name] !== nothing
    push!(obs_list, acc_obs[name])
    end
    if rv_obs[name] !== nothing
    push!(obs_list, rv_obs[name])
    end
    if get_data_flag(cfg, name, "z_prior") && zp_obs[name] !== nothing
        push!(obs_list, zp_obs[name])
    end
    if get_data_flag(cfg, name, "escape_velocity") && ev_obs[name] !== nothing
        push!(obs_list, ev_obs[name])
    end

    star = Planet(
        name = name,
        basis = Visual{KepOrbit},
        observations = obs_list,
        variables = @variables begin
            M = system.M
            P ~ P_prior                  # Period [yrs]
            a = cbrt(M * P^2)            # Semi-major axis [AU]
            e ~ e_prior                  # Eccentricity
            i ~ i_prior                  # Inclination [rad]
            ω ~ ω_prior                  # Argument of periastron [rad]
            Ω ~ Ω_prior                  # Longitude of ascending node [rad]
            θ ~ θ_prior                  # Mean anomaly at reference epoch [rad]
            tp = θ_at_epoch_to_tperi(θ, $epoch_mjd; a=a, e=e, i=i, ω=ω, Ω=Ω, M=M)
        end
    )
    push!(companions, star)
end

# === 4. Define the full system ===
sys_priors = cfg["priors"]["system"]
plx_prior     = parse_prior(sys_priors["plx"])
M_prior       = parse_prior(sys_priors["M"])
offsetx_prior = parse_prior(sys_priors["offsetx"])
offsety_prior = parse_prior(sys_priors["offsety"])

sys = System(
    name = get(cfg["meta"], "system_name", "Omega_Cen"),
    observations = [],
    companions = companions,
    variables = @variables begin
        plx ~ plx_prior              # Parallax [mas]
        M ~ M_prior                  # Host mass [solar masses]
        offsetx ~ offsetx_prior      # IMBH RA offset from assumed center [mas]
        offsety ~ offsety_prior      # IMBH Dec offset from assumed center [mas]
    end
)

# === 5. Model ===
model = Octofitter.LogDensityModel(sys)

# === 6. Sampling config ===
sampling_cfg         = cfg["sampling"]
n_rounds             = sampling_cfg["n_rounds"]
n_chains             = sampling_cfg["n_chains"]
n_chains_variational = sampling_cfg["n_chains_variational"]
checkpoint           = get(sampling_cfg, "checkpoint", false)

# === 7. Output config ===
slurm_job_id = get(ENV, "SLURM_JOB_ID", "none")
paths_cfg    = cfg["paths"]
output_dir   = isabspath(paths_cfg["output_dir"]) ? paths_cfg["output_dir"] : joinpath(dirname(abspath(config_path)), paths_cfg["output_dir"])
stars_tag    = join(star_names, "")
run_prefix   = "stars$(stars_tag)_$(n_chains)c_$(n_rounds)r$(is_resume ? "_cont" : "")_$(slurm_job_id)"
mkpath(output_dir)

# === 8. Write run summary ===
summary_path = joinpath(output_dir, "$(run_prefix)_summary.md")
open(summary_path, "w") do io
    println(io, "# Run Summary")
    println(io)
    println(io, "- **Date:** $(Dates.now())")
    println(io, "- **Slurm Job ID:** $(slurm_job_id)")
    println(io, "- **Stars:** $(join(star_names, ", "))")
    println(io, "- **Reference epoch:** $(epoch_mjd) MJD ($(epoch_year) yr)")
    println(io, "- **Config file:** $(abspath(config_path))")
    if is_resume
        println(io, "- **Run type:** continuation (target: $(n_rounds) total rounds)")
        println(io, "- **Resumed from PT checkpoint:** $(abspath(resume_pt_folder))")
    else
        println(io, "- **Run type:** fresh")
    end
    println(io)
    println(io, "## Sampling Parameters")
    println(io)
    println(io, "| Parameter | Value |")
    println(io, "|---|---|")
    println(io, "| n_rounds | $(n_rounds)$(is_resume ? " (total)" : "") |")
    println(io, "| n_chains | $(n_chains) |")
    println(io, "| n_chains_variational | $(n_chains_variational) |")
    println(io, "| checkpoint | $(checkpoint) |")
    println(io)
    println(io, "## System Priors")
    println(io)
    println(io, "| Parameter | Prior |")
    println(io, "|---|---|")
    for (k, v) in sys_priors
        println(io, "| $(k) | $(v) |")
    end
    if z_prior_sigma !== nothing
        println(io, "| z_prior | Normal(0, $(z_prior_sigma)) AU |")
    end
    println(io)
    println(io, "## Companion Priors (defaults)")
    println(io)
    println(io, "| Parameter | Prior |")
    println(io, "|---|---|")
    for (k, v) in cfg["priors"]["companion_defaults"]
        println(io, "| $(k) | $(v) |")
    end
    # Show per-star overrides if any
    overrides = get(get(cfg, "priors", Dict()), "overrides", Dict())
    if !isempty(overrides)
        println(io)
        println(io, "## Per-Star Prior Overrides")
        println(io)
        for (star, params) in overrides
            for (k, v) in params
                println(io, "- **Star $(star)**: $(k) = $(v)")
            end
        end
    end
    println(io)
    println(io, "## Full Configuration")
    println(io)
    println(io, "```toml")
    println(io, read(config_path, String))
    println(io, "```")
end
println("Run summary written to $(summary_path)")

# === 9. Fit with Pigeons ===
n_additional = nothing   # set below when resuming; used in summary
if is_resume
    println("Loading PT checkpoint from: $resume_pt_folder")
    pt_prev        = Pigeons.PT(resume_pt_folder)
    n_rounds_done  = pt_prev.inputs.n_rounds
    n_additional   = n_rounds - n_rounds_done
    n_additional > 0 || error(
        "n_rounds=$(n_rounds) ≤ rounds already completed ($(n_rounds_done)). " *
        "Increase n_rounds in [sampling] to a value greater than $(n_rounds_done)."
    )
    println("Resuming: $(n_rounds_done) rounds done, running $(n_additional) additional → $(n_rounds) total")
    Pigeons.increment_n_rounds!(pt_prev, n_additional)
    chain_pt = octofit_pigeons(pt_prev)
else
    chain_pt = octofit_pigeons(model; n_rounds, n_chains, n_chains_variational, checkpoint)
end
chain = chain_pt.chain
pt    = chain_pt.pt
println(chain)

# Record the PT exec folder so checkpoint files can be located and deleted later.
pt_exec_folder = try; pt.exec_folder; catch _
    try; pt.shared.exec_folder; catch _; "unknown"; end
end
println("PT exec folder (intermediate checkpoints): $pt_exec_folder")

# Write a small file containing just the PT exec folder path.
# A subsequent run can locate this by job_id to find the checkpoint for resuming.
pt_location_path = joinpath(output_dir, "$(run_prefix)_pt_location.txt")
write(pt_location_path, pt_exec_folder)
println("PT location written to: $pt_location_path")

open(summary_path, "a") do io
    println(io)
    println(io, "## Sampling Result")
    println(io)
    if is_resume && n_additional !== nothing
        println(io, "- **Additional rounds run:** $(n_additional)")
        println(io, "- **Total rounds:** $(n_rounds)")
    end
    println(io)
    println(io, "## Pigeons PT Checkpoint")
    println(io)
    println(io, "The PT exec folder contains intermediate checkpoint files.")
    println(io, "It can be deleted once the chain FITS file has been verified.")
    println(io)
    println(io, "| | Path |")
    println(io, "|---|---|")
    println(io, "| PT exec folder | `$(pt_exec_folder)` |")
    println(io, "| PT location file | `$(pt_location_path)` |")
end

# === 10. Save Chain ===
Octofitter.savechain(joinpath(output_dir, "$(run_prefix)_chain.fits"), chain)

# === 11. Generate plots ===
ARGS_bak = copy(ARGS)
empty!(ARGS)
push!(ARGS, joinpath(output_dir, "$(run_prefix)_chain.fits"))
include(joinpath(@__DIR__, "plot_chain.jl"))
empty!(ARGS)
append!(ARGS, ARGS_bak)
