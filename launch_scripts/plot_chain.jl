"""
Generate plots from a saved chain FITS file.

Re-creates the corner plot, sky-plane orbit panels, and posterior histogram
panels that are normally produced at the end of octo_orbit_direct_likelihoods.jl.

Usage (standalone, safe on a login node):
    julia --project=<Octofitter_imbh.jl> plot_chain.jl <chain.fits>

When called via include() from octo_orbit_direct_likelihoods.jl, the thread
count of the parent process is inherited automatically.

The companion *_summary.md file (same directory, same run prefix) is loaded
automatically to recover the epoch, priors, and star list.  The full TOML
configuration is embedded in that file, so no separate config.toml is needed.

Star names are also cross-checked against the chain column names (columns
ending in _a with matching _e, _i, _ω, _Ω, _tp companions).

Output files are written to the same directory as *chain.fits, using the same
run_prefix (everything before _chain.fits in the filename).
"""

ENV["OCTOFITTERPY_AUTOLOAD_EXTENSIONS"] = "yes"

using Octofitter
using Octofitter: @variables, System
using CairoMakie
using PairPlots
using Distributions
using Unitful
using UnitfulAstro
using LinearAlgebra
using Statistics
using Printf
using TOML

push!(LOAD_PATH, @__DIR__)
using octo_utils

include(joinpath(@__DIR__, "parse_config.jl"))

# ── 1. Parse arguments ───────────────────────────────────────────────────────

length(ARGS) >= 1 || error("Usage: julia plot_chain.jl <chain.fits>")
chain_path = ARGS[1]
isfile(chain_path) || error("Chain file not found: $chain_path")

output_dir = dirname(abspath(chain_path))
chain_basename = basename(chain_path)
run_prefix = endswith(chain_basename, "_chain.fits") ?
    chain_basename[1:end-length("_chain.fits")] : splitext(chain_basename)[1]

# ── 2. Load configuration from the companion summary.md ──────────────────────

summary_path = joinpath(output_dir, "$(run_prefix)_summary.md")
isfile(summary_path) || error(
    "Summary file not found: $summary_path\n" *
    "Expected a *_summary.md alongside the chain file.")

summary_text = read(summary_path, String)

# Extract the TOML block embedded between ```toml and ``` fences
toml_match = match(r"```toml\r?\n(.*?)```"s, summary_text)
toml_match === nothing && error(
    "Could not find a ```toml ... ``` block in $summary_path")
cfg = TOML.parse(toml_match[1])
println("Loaded configuration from: $summary_path")

epoch_mjd  = get_epoch_mjd(cfg)
epoch_year = cfg["epoch"]["year"]
println("Reference epoch: $epoch_mjd MJD ($epoch_year yr)")

# ── 3. Load chain ────────────────────────────────────────────────────────────

chain = Octofitter.loadchain(chain_path)
println("Loaded chain: $chain_path")
println(chain)

# ── 4. Get star names from summary, verify against chain columns ──────────────
# FITS column names are ASCII-only; Unicode characters (ω, Ω) may be encoded
# differently on save/load.  Use the summary as the authoritative source for
# star names, then discover the actual column names for each orbital element.

col_names = Set(Symbol.(names(chain)))

# Primary: parse from summary
summary_stars_line = match(r"\*\*Stars:\*\*\s*([^\n]+)", summary_text)
summary_stars_line !== nothing ||
    error("Could not find '**Stars:**' line in $summary_path")
star_names = sort!(String.(strip.(split(summary_stars_line[1], ","))))
println("Star names from summary: $(join(star_names, ", "))")

# Helper: find a chain column for a given star and orbital element, trying
# multiple name variants to handle FITS ASCII encoding of Unicode symbols.
const _col_variants = Dict(
    "a"  => ["a"],
    "e"  => ["e"],
    "i"  => ["i"],
    "ω"  => ["ω",  "omega", "w"],
    "Ω"  => ["Ω",  "Omega", "W"],
    "tp" => ["tp"],
)
function find_col(col_names, star, element)
    for v in _col_variants[element]
        sym = Symbol("$(star)_$(v)")
        sym in col_names && return sym
    end
    error("Cannot find chain column for $(star)_$(element). " *
          "Available columns with prefix '$(star)_': " *
          join(filter(c -> startswith(String(c), "$(star)_"), collect(col_names)), ", "))
end

# Verify all expected columns exist
for name in star_names, el in keys(_col_variants)
    find_col(col_names, name, el)   # errors early with a clear message if missing
end

# ── 5. Rebuild observation objects ───────────────────────────────────────────

astrom_obs = Dict{String, Any}()
pm_obs     = Dict{String, Any}()
acc_obs    = Dict{String, Any}()
rv_obs     = Dict{String, Any}()   # may hold nothing for stars without RV
for name in star_names
    haskey(octo_utils.stars, name) ||
        error("Star '$name' not found in octo_utils.stars.")
    star = octo_utils.stars[name]
    a, p, ac, rv, _zp = octo_utils.build_star_observations(star, epoch_mjd)
    astrom_obs[name] = a
    pm_obs[name]     = p
    acc_obs[name]    = ac
    rv_obs[name]     = rv
end

# ── 6. Rebuild the Octofitter model (required by octocorner) ─────────────────

companions = Planet[]
for name in star_names
    P_prior = parse_prior(get_companion_prior(cfg, name, "P"))
    e_prior = parse_prior(get_companion_prior(cfg, name, "e"))
    i_prior = parse_prior(get_companion_prior(cfg, name, "i"))
    ω_prior = parse_prior(get_companion_prior(cfg, name, "omega"))
    Ω_prior = parse_prior(get_companion_prior(cfg, name, "Omega"))
    θ_prior = parse_prior(get_companion_prior(cfg, name, "theta"))

    star = Planet(
        name = name,
        basis = Visual{KepOrbit},
        observations = [ObsPriorAstromONeil2019(astrom_obs[name]), pm_obs[name], acc_obs[name]],
        variables = @variables begin
            M = system.M
            P ~ P_prior
            a = cbrt(M * P^2)
            e ~ e_prior
            i ~ i_prior
            ω ~ ω_prior
            Ω ~ Ω_prior
            θ ~ θ_prior
            tp = θ_at_epoch_to_tperi(θ, $epoch_mjd; a=a, e=e, i=i, ω=ω, Ω=Ω, M=M)
        end
    )
    push!(companions, star)
end

sys_priors    = cfg["priors"]["system"]
plx_prior     = parse_prior(sys_priors["plx"])
M_prior       = parse_prior(sys_priors["M"])
offsetx_prior = parse_prior(sys_priors["offsetx"])
offsety_prior = parse_prior(sys_priors["offsety"])

sys = System(
    name = get(cfg["meta"], "system_name", "Omega_Cen"),
    observations = [],
    companions = companions,
    variables = @variables begin
        plx ~ plx_prior
        M ~ M_prior
        offsetx ~ offsetx_prior
        offsety ~ offsety_prior
    end
)

model = Octofitter.LogDensityModel(sys)

# ── 7. Extract posterior samples ─────────────────────────────────────────────

M_samples   = vec(chain[:M])
plx_samples = vec(chain[:plx])
ox_samples  = vec(chain[:offsetx])
oy_samples  = vec(chain[:offsety])

star_samples = Dict{String, NamedTuple}()
for name in star_names
    star_samples[name] = (
        a  = vec(chain[find_col(col_names, name, "a")]),
        e  = vec(chain[find_col(col_names, name, "e")]),
        i  = vec(chain[find_col(col_names, name, "i")]),
        ω  = vec(chain[find_col(col_names, name, "ω")]),
        Ω  = vec(chain[find_col(col_names, name, "Ω")]),
        tp = vec(chain[find_col(col_names, name, "tp")]),
    )
end

# ── 8. Posterior summaries ───────────────────────────────────────────────────

function format_stat(label, samples; scale=1.0)
    med = median(samples) * scale
    lo  = quantile(samples, 0.16) * scale
    hi  = quantile(samples, 0.84) * scale
    return @sprintf("%-20s  %10.3f  [%8.3f, %8.3f]", label, med, lo, hi)
end

stat_lines = String[]
push!(stat_lines, @sprintf("%-20s  %10s  [%8s, %8s]", "Param", "Median", "16%", "84%"))
push!(stat_lines, format_stat("M_IMBH [10⁴ M☉]", M_samples; scale=1e-4))
push!(stat_lines, format_stat("plx [mas]",        plx_samples))
push!(stat_lines, format_stat("offsetx [mas]",    ox_samples))
push!(stat_lines, format_stat("offsety [mas]",    oy_samples))
for name in star_names
    s = star_samples[name]
    push!(stat_lines, format_stat("$(name): a [AU]", s.a))
    push!(stat_lines, format_stat("$(name): e",      s.e))
    push!(stat_lines, format_stat("$(name): i [°]",  s.i; scale=180/π))
    push!(stat_lines, format_stat("$(name): ω [°]",  s.ω; scale=180/π))
    push!(stat_lines, format_stat("$(name): Ω [°]",  s.Ω; scale=180/π))
    push!(stat_lines, format_stat("$(name): tp [mjd]", s.tp;))
end

println("\n=== Posterior summaries (median, 68% CI) ===")
for line in stat_lines
    println(line)
end
# Note: _posterior_stats.txt is written at the very end of the script so that
# physical-plausibility diagnostics (Sections 11.5–11.7) can append to
# stat_lines before the file is produced.
stats_path = joinpath(output_dir, "$(run_prefix)_posterior_stats.txt")

# ── 9. Corner plot ───────────────────────────────────────────────────────────

println("\nGenerating corner plot...")
# Subsample chain for corner plot to limit memory usage
max_corner_samples = 10_000
if size(chain, 1) > max_corner_samples
    thin_idx = round.(Int, range(1, size(chain, 1), length=max_corner_samples))
    chain_thin = chain[thin_idx, :, :]
else
    chain_thin = chain
end
corner_plot = octocorner(model, chain_thin; small=true,
    includecols=["M", "offsetx", "offsety"],
    labels=Dict{Symbol,Any}(
        :offsetx => "Δα*_IMBH [mas]",
        :offsety => "Δδ_IMBH [mas]",
    )
)
save(joinpath(output_dir, "$(run_prefix)_corner.png"), corner_plot, px_per_unit=3)
corner_plot = nothing; chain_thin = nothing; GC.gc()
println("Corner plot saved.")

# ── 10. Sky-plane orbit panels ───────────────────────────────────────────────

println("Generating orbit panels...")

# Per-star color: cycle through Wong colors, skipping the first (blue → too similar to black).
star_colors = Dict{String, Any}(
    name => Makie.wong_colors()[mod1(k + 1, length(Makie.wong_colors()))]
    for (k, name) in enumerate(star_names)
)

# Fraction of auto-computed range to display per star, centred on IMBH median
const star_zoom = Dict("A" => 0.60, "C" => 0.80, "D" => 0.50, "E" => 0.50, "F" => 0.50)

sample_idx = round.(Int, range(1, length(M_samples), length=100))

function star_orbit_panel!(ax, s, M_samp, plx_samp, ox_samp, oy_samp,
                            obs_ra, obs_dec, obs_pmra, obs_pmdec,
                            epoch_mjd, sample_idx, color;
                            scale_pm=250.0)
    ox_med_loc = median(ox_samp)
    oy_med_loc = median(oy_samp)
    for idx in sample_idx
        orb_s = Visual{KepOrbit}(;
            a=s.a[idx], e=s.e[idx], i=s.i[idx],
            ω=s.ω[idx], Ω=s.Ω[idx], tp=s.tp[idx],
            M=M_samp[idx], plx=plx_samp[idx])
        P_s = s.a[idx]^1.5 / sqrt(M_samp[idx])  # period in years for this sample
        ts  = range(epoch_mjd, epoch_mjd + P_s * 365.25; length=300)
        ra_s  = [raoff(orbitsolve(orb_s, t)) + ox_samp[idx] for t in ts]
        dec_s = [decoff(orbitsolve(orb_s, t)) + oy_samp[idx] for t in ts]
        lines!(ax, ra_s, dec_s; color=(color, 0.5), linewidth=0.5)
    end
    scatter!(ax, [0.0], [0.0]; marker='+', markersize=20, color=:black)
    scatter!(ax, [ox_med_loc], [oy_med_loc]; marker=:circle, markersize=12, color=:black)
    
    
    arrows2d!(ax, [obs_ra], [obs_dec],
        [obs_pmra * scale_pm], [obs_pmdec * scale_pm];
        color=:red, shaftwidth=2.0, tipwidth=10, tiplength=10)
    scatter!(ax, [obs_ra], [obs_dec];
        marker='★', color=Makie.wong_colors()[2], markersize=14,
        strokecolor=:black, strokewidth=0.5)
end

n_stars    = length(star_names)
n_cols_orb = min(n_stars, 3)
n_rows_orb = ceil(Int, n_stars / n_cols_orb)
# If the last row of individual panels has empty cells, reuse them for the
# combined panel; otherwise place it on a new row.
n_filled_last  = mod(n_stars, n_cols_orb)   # 0 means last row is full
combined_row   = n_filled_last == 0 ? n_rows_orb + 1 : n_rows_orb
combined_cols  = n_filled_last == 0 ? (1:n_cols_orb) : ((n_filled_last + 1):n_cols_orb)
n_rows_fig     = n_filled_last == 0 ? n_rows_orb + 1 : n_rows_orb

fig_orbits = Figure(size=(n_cols_orb * 420, n_rows_fig * 420), fontsize=18)

for (k, name) in enumerate(star_names)
    row   = ceil(Int, k / n_cols_orb)
    col   = mod1(k, n_cols_orb)
    color = star_colors[name]
    ax    = Axis(fig_orbits[row, col];
        xlabel="Δα* [mas]", ylabel="Δδ [mas]",
        xreversed=true, aspect=DataAspect(),
        xgridvisible=false, ygridvisible=false)
    star_orbit_panel!(ax, star_samples[name], M_samples, plx_samples,
        ox_samples, oy_samples,
        astrom_obs[name].table.ra[1], astrom_obs[name].table.dec[1],
        pm_obs[name].table.pmra[1], pm_obs[name].table.pmdec[1],
        epoch_mjd, sample_idx, color)
    text!(ax, "Star $name"; position=(0.05, 0.95), align=(:left, :top),
          space=:relative, fontsize=18)
    # Zoom and centre on IMBH median position
    autolimits!(ax)
    fl   = ax.finallimits[]
    half = max(fl.widths[1], fl.widths[2]) / 2 * get(star_zoom, name, 1.0)
    ox_med = median(ox_samples)
    oy_med = median(oy_samples)
    limits!(ax, ox_med - half, ox_med + half, oy_med - half, oy_med + half)
    ax.xreversed[] = true   # re-assert after limits!() which can reset the direction
end

# ── Combined panel: all stars on one plate ───────────────────────────────────
ax_all = Axis(fig_orbits[combined_row, combined_cols];
    xlabel="Δα* [mas]", ylabel="Δδ [mas]",
    xreversed=true, aspect=DataAspect(),
    xgridvisible=false, ygridvisible=false)

for (k, name) in enumerate(star_names)
    color    = star_colors[name]
    s        = star_samples[name]
    obs_ra   = astrom_obs[name].table.ra[1]
    obs_dec  = astrom_obs[name].table.dec[1]
    for idx in sample_idx
        orb_s = Visual{KepOrbit}(;
            a=s.a[idx], e=s.e[idx], i=s.i[idx],
            ω=s.ω[idx], Ω=s.Ω[idx], tp=s.tp[idx],
            M=M_samples[idx], plx=plx_samples[idx])
        P_s = s.a[idx]^1.5 / sqrt(M_samples[idx])  # period in years for this sample
        ts  = range(epoch_mjd, epoch_mjd + P_s * 365.25; length=300)
        ra_s  = [raoff(orbitsolve(orb_s, t)) + ox_samples[idx] for t in ts]
        dec_s = [decoff(orbitsolve(orb_s, t)) + oy_samples[idx] for t in ts]
        lines!(ax_all, ra_s, dec_s; color=(color, 0.3), linewidth=0.5)
    end
    scatter!(ax_all, [obs_ra], [obs_dec];
        marker='★', color=Makie.wong_colors()[2], markersize=14,
        strokecolor=:black, strokewidth=0.5)
end

scatter!(ax_all, [0.0], [0.0]; marker='+', markersize=20, color=:black, label="AvdM10 centre")
scatter!(ax_all, [median(ox_samples)], [median(oy_samples)];
    marker=:circle, markersize=12, color=:black, label="IMBH")
axislegend(ax_all; position=:rt, framevisible=false)
autolimits!(ax_all)
fl_all = ax_all.finallimits[]
half_all = max(fl_all.widths[1], fl_all.widths[2]) / 2 * get(star_zoom, "F", 0.50)
limits!(ax_all, median(ox_samples) - half_all, median(ox_samples) + half_all,
                median(oy_samples) - half_all, median(oy_samples) + half_all)
ax_all.xreversed[] = true   # re-assert after limits!() which can reset the direction

save(joinpath(output_dir, "$(run_prefix)_orbit_panels.png"), fig_orbits, px_per_unit=3)
fig_orbits = nothing; GC.gc()
println("Orbit panels saved.")

# ── 11. Posterior histogram panels ───────────────────────────────────────────

println("Generating posterior panels...")

function param_panel!(layout, row, col, color, samples, xlabel;
                       show_legend=false, xlims=nothing, bins=30, xticks=Makie.automatic)
    ax = Axis(layout[row, col]; xlabel=xlabel, ylabel="Probability Density",
              xgridvisible=false, ygridvisible=false, xticks=xticks)
    med = median(samples)
    hist!(ax, samples; normalization=:pdf, bins=bins, color=(color, 0.7))
    vlines!(ax, [med]; color=Makie.wong_colors()[2], linestyle=:solid, label="Median")
    show_legend && axislegend(ax; position=:rt, framevisible=false)
    xlims !== nothing && Makie.xlims!(ax, xlims...)
end

fig_post = Figure(size=(1600, (1 + n_stars) * 260), fontsize=18)

sys_color = Makie.wong_colors()[1]
param_panel!(fig_post, 1, 1, sys_color, M_samples ./ 1e4,
    Makie.rich("M", Makie.subscript("IMBH"), " [10⁴ M", Makie.subscript("☉"), "]");
    show_legend=true)
param_panel!(fig_post, 1, 2, sys_color, plx_samples, "plx [mas]")
param_panel!(fig_post, 1, 3, sys_color, ox_samples,
    Makie.rich("Δα*", Makie.subscript("IMBH"), " [mas]"))
param_panel!(fig_post, 1, 4, sys_color, oy_samples,
    Makie.rich("Δδ", Makie.subscript("IMBH"), " [mas]"))

for (k, name) in enumerate(star_names)
    row  = k + 1
    s    = star_samples[name]
    c    = star_colors[name]
    param_panel!(fig_post, row, 1, c, s.a,           "$(name): a [AU]";
        xlims=(0, 20_000),
        bins=range(0, 20_000, length=51),
        xticks=[0, 10_000, 20_000])
    param_panel!(fig_post, row, 2, c, s.e,           "$(name): e")
    param_panel!(fig_post, row, 3, c, rad2deg.(s.i), "$(name): i [°]")
    param_panel!(fig_post, row, 4, c, rad2deg.(s.ω), "$(name): ω [°]")
    param_panel!(fig_post, row, 5, c, rad2deg.(s.Ω), "$(name): Ω [°]")
end
save(joinpath(output_dir, "$(run_prefix)_posteriors.png"), fig_post, px_per_unit=3)
fig_post = nothing; GC.gc()
println("Posterior panels saved.")

# ── 11.5. Physical plausibility: pericenter, velocity, tidal radii ───────────

println("Generating plausibility diagnostics...")

# 1 AU/yr in km/s (1 AU / 1 yr in SI)
const AU_YR_TO_KMS = 4.7404705
# G·M_sun in AU³/yr² (Kepler's third law with a in AU, P in yr, M in M_sun)
const FOUR_PI2 = 4 * π^2

# Reference stellar templates (not priors — for tidal radius lines only)
const R_SUN_AU   = 1 / 215.032           # 1 R_sun in AU
const R_GIANT_AU = 30 * R_SUN_AU         # ~30 R_sun for a red giant
tidal_radius(R_star_au, m_star_Msun, M_BH_Msun) =
    R_star_au * cbrt(M_BH_Msun / m_star_Msun)

"""
Per-sample orbital scalars for one star: pericenter/apocenter distances (AU),
pericenter/apocenter speeds (km/s, via vis-viva), and period (yr).
"""
function orbital_scalars(s, M_samples)
    a = s.a; e = s.e
    r_peri = a .* (1 .- e)
    r_apo  = a .* (1 .+ e)
    v_peri = @. sqrt(FOUR_PI2 * M_samples * (1 + e) / (a * (1 - e))) * AU_YR_TO_KMS
    v_apo  = @. sqrt(FOUR_PI2 * M_samples * (1 - e) / (a * (1 + e))) * AU_YR_TO_KMS
    P_yr   = @. sqrt(a^3 / M_samples)
    return (; r_peri, r_apo, v_peri, v_apo, P_yr)
end

# Mass-dependent reference radii (one value per posterior draw)
rt_ms_samples  = tidal_radius.(R_SUN_AU,   1.0, M_samples)
rt_rg_samples  = tidal_radius.(R_GIANT_AU, 0.8, M_samples)
# Schwarzschild radius in AU:  r_s = 2GM/c² ≈ M[M_sun] · 1.909e-8 AU
r_schw_samples = M_samples .* 1.909e-8

# Append scalars to the text summary
for name in star_names
    scal = orbital_scalars(star_samples[name], M_samples)
    push!(stat_lines, format_stat("$(name): r_peri [AU]", scal.r_peri))
    push!(stat_lines, format_stat("$(name): r_apo  [AU]", scal.r_apo))
    push!(stat_lines, format_stat("$(name): v_peri [km/s]", scal.v_peri))
    push!(stat_lines, format_stat("$(name): v_apo  [km/s]", scal.v_apo))
    push!(stat_lines, format_stat("$(name): P      [yr]",   scal.P_yr))
end
push!(stat_lines, format_stat("r_tidal MS  [AU]", rt_ms_samples))
push!(stat_lines, format_stat("r_tidal RG  [AU]", rt_rg_samples))
push!(stat_lines, format_stat("r_Schw      [AU]", r_schw_samples))

fig_phys = Figure(size=(1200, n_stars * 260), fontsize=18)
for (k, name) in enumerate(star_names)
    s    = star_samples[name]
    scal = orbital_scalars(s, M_samples)
    c    = star_colors[name]

    ax_r = Axis(fig_phys[k, 1];
        xlabel="$(name): r_peri [AU]", ylabel="Probability Density",
        xscale=log10,
        xgridvisible=false, ygridvisible=false)
    # Build log-spaced bins so the histogram renders correctly on a log axis
    r_lo = max(minimum(scal.r_peri),
               0.5 * min(median(rt_ms_samples), median(r_schw_samples)))
    r_hi = maximum(scal.r_peri)
    r_bins = 10 .^ range(log10(r_lo), log10(r_hi); length=31)
    hist!(ax_r, scal.r_peri; normalization=:pdf, bins=r_bins, color=(c, 0.7))
    vlines!(ax_r, [median(rt_ms_samples)];
            color=:steelblue, linestyle=:dash, label="r_t (MS)")
    vlines!(ax_r, [median(rt_rg_samples)];
            color=:firebrick, linestyle=:dash, label="r_t (RG)")
    vlines!(ax_r, [median(r_schw_samples)];
            color=:black, linestyle=:dot, label="r_Schw")
    k == 1 && axislegend(ax_r; position=:lt, framevisible=false)

    ax_v = Axis(fig_phys[k, 2];
        xlabel="$(name): v_peri [km/s]", ylabel="Probability Density",
        xgridvisible=false, ygridvisible=false)
    hist!(ax_v, scal.v_peri; normalization=:pdf, bins=30, color=(c, 0.7))
end
save(joinpath(output_dir, "$(run_prefix)_plausibility.png"), fig_phys, px_per_unit=3)
fig_phys = nothing; GC.gc()
println("Plausibility diagnostics saved.")

# ── 11.6. True anomaly at obs epoch + acceleration-vector alignment ──────────

println("Generating phase / acceleration-alignment diagnostics...")

using PlanetOrbits: trueanom, radvel, accra, accdec

"True anomaly (deg, wrapped to (-180, 180]) at epoch_mjd, one per chain draw."
function true_anomaly_at_epoch(s, M_samp, plx_samp, epoch_mjd)
    ν = Vector{Float64}(undef, length(M_samp))
    @inbounds for idx in eachindex(M_samp)
        orb = Visual{KepOrbit}(;
            a=s.a[idx], e=s.e[idx], i=s.i[idx],
            ω=s.ω[idx], Ω=s.Ω[idx], tp=s.tp[idx],
            M=M_samp[idx], plx=plx_samp[idx])
        ν[idx] = rad2deg(trueanom(orbitsolve(orb, epoch_mjd)))
    end
    return @. mod(ν + 180, 360) - 180
end

"Angle (deg, 0–180) between measured accel vector and star→IMBH direction."
function accel_alignment_angle(name, ox_samp, oy_samp)
    obs_ra  = astrom_obs[name].table.ra[1]
    obs_dec = astrom_obs[name].table.dec[1]
    ax_meas = acc_obs[name].table.accra[1]
    ay_meas = acc_obs[name].table.accdec[1]
    a_norm  = hypot(ax_meas, ay_meas)
    ax_hat  = ax_meas / a_norm
    ay_hat  = ay_meas / a_norm
    dx = ox_samp .- obs_ra
    dy = oy_samp .- obs_dec
    r  = hypot.(dx, dy)
    dxh = dx ./ r
    dyh = dy ./ r
    cosφ = clamp.(ax_hat .* dxh .+ ay_hat .* dyh, -1.0, 1.0)
    return rad2deg.(acos.(cosφ))
end

"""
    accel_toward_imbh_zscore(name, ox_samp, oy_samp)

Per-posterior-sample z-score: the component of the measured acceleration in
the direction of the star→IMBH unit vector, divided by the propagated
uncertainty of that component.  z > 0 means the measured acceleration has a
component pointing toward the IMBH; z ≈ 0 means the measurement cannot
distinguish the IMBH direction; z < 0 means it points away.

Also returns the angular uncertainty on the measured acceleration direction
(in degrees), which is constant across the posterior.
"""
function accel_toward_imbh_zscore(name, ox_samp, oy_samp)
    obs_ra  = astrom_obs[name].table.ra[1]
    obs_dec = astrom_obs[name].table.dec[1]
    ax_meas = acc_obs[name].table.accra[1]
    ay_meas = acc_obs[name].table.accdec[1]
    σx      = acc_obs[name].table.σ_accra[1]
    σy      = acc_obs[name].table.σ_accdec[1]

    # Angular uncertainty of the measured direction via error propagation on atan2
    σ_φ_deg = rad2deg(sqrt((ay_meas * σx)^2 + (ax_meas * σy)^2) / (ax_meas^2 + ay_meas^2))

    # Per-sample unit vector from star toward IMBH
    dx = ox_samp .- obs_ra
    dy = oy_samp .- obs_dec
    r  = hypot.(dx, dy)
    dxh = dx ./ r
    dyh = dy ./ r

    # Component of measured acceleration in the IMBH direction and its uncertainty
    a_toward  = ax_meas .* dxh .+ ay_meas .* dyh
    σ_toward  = sqrt.((σx .* dxh).^2 .+ (σy .* dyh).^2)

    z = a_toward ./ σ_toward
    return z, σ_φ_deg
end

fig_pa = Figure(size=(1500, n_stars * 240), fontsize=18)
for (k, name) in enumerate(star_names)
    c  = star_colors[name]
    ν  = true_anomaly_at_epoch(star_samples[name], M_samples, plx_samples, epoch_mjd)
    Δφ = accel_alignment_angle(name, ox_samples, oy_samples)
    z_accel, σ_φ_deg = accel_toward_imbh_zscore(name, ox_samples, oy_samples)

    ax_ν = Axis(fig_pa[k, 1];
        xlabel="$(name): ν(t_obs) [°]", ylabel="Probability Density",
        xticks=-180:90:180,
        xgridvisible=false, ygridvisible=false)
    hist!(ax_ν, ν; normalization=:pdf, bins=40, color=(c, 0.7))
    vlines!(ax_ν, [0.0]; color=:black, linestyle=:dot)

    ax_φ = Axis(fig_pa[k, 2];
        xlabel="$(name): accel misalignment Δφ [°]",
        ylabel="Probability Density",
        xgridvisible=false, ygridvisible=false)
    hist!(ax_φ, Δφ; normalization=:pdf, bins=40, color=(c, 0.7))
    # Dashed line at 0° (perfect alignment) and shaded band showing the
    # angular uncertainty of the measured acceleration direction
    vlines!(ax_φ, [0.0]; color=:black, linestyle=:dot, label="Perfect alignment")
    vspan!(ax_φ, 0.0, σ_φ_deg; color=(:grey, 0.25), label="Meas. dir. uncertainty (1σ)")
    k == 1 && axislegend(ax_φ; position=:rt, framevisible=false, labelsize=12)

    # Column 3: z-score (component of measured acc toward IMBH / uncertainty)
    ax_z = Axis(fig_pa[k, 3];
        xlabel="$(name): accel z-score toward IMBH",
        ylabel="Probability Density",
        xgridvisible=false, ygridvisible=false)
    hist!(ax_z, z_accel; normalization=:pdf, bins=40, color=(c, 0.7))
    vlines!(ax_z, [0.0]; color=:black, linestyle=:dot)
    vlines!(ax_z, [-1.0, 1.0]; color=:grey, linestyle=:dash)
    k == 1 && text!(ax_z, "z>0: acc toward IMBH"; position=(0.55, 0.90),
                    align=(:left, :top), space=:relative, fontsize=11, color=:grey40)

    push!(stat_lines, format_stat("$(name): ν(t_obs) [°]", ν))
    push!(stat_lines, format_stat("$(name): Δφ_accel [°]", Δφ))
    push!(stat_lines,
          @sprintf("%-20s  %10.1f °", "$(name): σ_φ_acc (meas)", σ_φ_deg))
    push!(stat_lines, format_stat("$(name): z_acc→IMBH", z_accel))
end
save(joinpath(output_dir, "$(run_prefix)_phase_accel.png"), fig_pa, px_per_unit=3)
fig_pa = nothing; GC.gc()
println("Phase / acceleration alignment diagnostics saved.")

# ── 11.7. Radial-velocity consistency check (stars with RV data only) ────────

rv_stars = [n for n in star_names if rv_obs[n] !== nothing]
if !isempty(rv_stars)
    println("Generating RV consistency check...")
    fig_rv = Figure(size=(500 * length(rv_stars), 400), fontsize=18)
    for (k, name) in enumerate(rv_stars)
        s = star_samples[name]
        rv_pred = Vector{Float64}(undef, length(M_samples))
        @inbounds for idx in eachindex(M_samples)
            orb = Visual{KepOrbit}(;
                a=s.a[idx], e=s.e[idx], i=s.i[idx],
                ω=s.ω[idx], Ω=s.Ω[idx], tp=s.tp[idx],
                M=M_samples[idx], plx=plx_samples[idx])
            rv_pred[idx] = radvel(orbitsolve(orb, epoch_mjd))  # m/s, peculiar
        end
        rv_pred_kms = rv_pred ./ 1000.0

        rv_meas_kms  = rv_obs[name].table.rv[1]   / 1000.0
        rv_sigma_kms = rv_obs[name].table.σ_rv[1] / 1000.0

        ax = Axis(fig_rv[1, k];
            xlabel="$(name): peculiar RV [km/s]",
            ylabel="Probability Density",
            xgridvisible=false, ygridvisible=false)
        hist!(ax, rv_pred_kms; normalization=:pdf, bins=40,
              color=(star_colors[name], 0.7), label="Posterior prediction")
        vspan!(ax, rv_meas_kms - rv_sigma_kms, rv_meas_kms + rv_sigma_kms;
               color=(:grey, 0.35), label="Measured ± 1σ")
        vlines!(ax, [rv_meas_kms]; color=:black, linewidth=2, label="Measured")
        k == 1 && axislegend(ax; position=:rt, framevisible=false)

        z = (median(rv_pred_kms) - rv_meas_kms) / rv_sigma_kms
        push!(stat_lines,
              @sprintf("%-20s  %+10.2f σ", "$(name): RV residual", z))
    end
    save(joinpath(output_dir, "$(run_prefix)_rv_check.png"), fig_rv, px_per_unit=3)
    fig_rv = nothing; GC.gc()
    println("RV consistency check saved.")
else
    println("No stars with RV data; skipping RV consistency check.")
end

# ── 11.8. Acceleration posterior predictive check ────────────────────────────
#
# Goal: treat the acceleration measurements as an independent validation of the
# orbit model rather than a fitting constraint.  For each posterior draw we
# compute the sky-plane acceleration that the fitted Keplerian orbit predicts
# at the observation epoch, then compare the full predictive distribution to the
# measured values.
#
# This is a standard posterior predictive check (PPC): if the model is
# consistent with the acceleration data, the measured acceleration should fall
# in the bulk of the predictive cloud.  If it sits in a clear tail, the
# acceleration is in tension with the model — which may indicate either a
# genuine physical inconsistency or a systematic error in the measurement.
#
# Two complementary diagnostics per star:
#
#   Left panel — 2D predictive scatter in (accra, accdec) space [mas/yr²]:
#     Each point is the sky-plane acceleration predicted by one posterior draw,
#     computed via accra(sol) and accdec(sol) from PlanetOrbits.  The measured
#     value is shown as a black cross with ±1σ error bars.  If the cross falls
#     in the bulk of the scatter cloud the model is consistent with the
#     acceleration data; if it is a clear outlier the two are in tension.
#
#   Right panel — Distribution of 2D chi-squared residuals:
#     For each posterior draw i the 2D chi-squared distance between the
#     predicted and measured acceleration is:
#
#         χ²_i = (accra_meas  - accra_pred_i )² / σ_ra²
#               + (accdec_meas - accdec_pred_i)² / σ_dec²
#
#     The vertical dashed line marks χ² = 2.30, which is the 68th percentile
#     of a chi-squared distribution with 2 degrees of freedom — i.e. the
#     boundary of the 1σ error ellipse in 2D.  The fraction of posterior draws
#     with χ²_i ≤ 2.30 is labelled f₆₈ on the plot:
#       • f₆₈ ≈ 0.68 → the measurement lies in the typical 1σ bulk of the
#         predictive distribution (fully consistent)
#       • f₆₈ << 0.68 → the measurement is in the tail; the model cannot
#         easily reproduce the observed acceleration, indicating tension
#
# Note: accra(sol) and accdec(sol) from PlanetOrbits return the gravitational
# acceleration of the star toward the IMBH in the sky plane (mas/yr²), with
# the same sign convention as PlanetAccelObs.  The IMBH offset (offsetx,
# offsety) is already encoded in the orbital geometry, so no additional
# correction is needed.

println("Generating acceleration posterior predictive check...")

fig_acc_ppc = Figure(size=(1000, n_stars * 280), fontsize=18)

for (k, name) in enumerate(star_names)
    s    = star_samples[name]
    c    = star_colors[name]

    # Measured acceleration and uncertainties for this star [mas/yr²]
    ax_meas = acc_obs[name].table.accra[1]
    ay_meas = acc_obs[name].table.accdec[1]
    σ_ra    = acc_obs[name].table.σ_accra[1]
    σ_dec   = acc_obs[name].table.σ_accdec[1]

    # --- Compute predicted acceleration for every posterior draw ---
    # accra(sol) and accdec(sol) evaluate the Keplerian gravitational
    # acceleration at the solved orbital position (mas/yr²).  This is the
    # model quantity that PlanetAccelObs compares to the measurement in the
    # likelihood — here we compute it purely as a prediction check.
    n_samp      = length(M_samples)
    accra_pred  = Vector{Float64}(undef, n_samp)
    accdec_pred = Vector{Float64}(undef, n_samp)
    @inbounds for idx in 1:n_samp
        orb = Visual{KepOrbit}(;
            a=s.a[idx], e=s.e[idx], i=s.i[idx],
            ω=s.ω[idx], Ω=s.Ω[idx], tp=s.tp[idx],
            M=M_samples[idx], plx=plx_samples[idx])
        sol = orbitsolve(orb, epoch_mjd)
        accra_pred[idx]  = accra(sol)
        accdec_pred[idx] = accdec(sol)
    end

    # --- 2D chi-squared residual per draw ---
    # Measures the distance (in units of measurement uncertainty) between
    # each predicted acceleration and the measured value.
    chi2 = @. (ax_meas - accra_pred)^2 / σ_ra^2 +
               (ay_meas - accdec_pred)^2 / σ_dec^2

    # Fraction of draws that predict an acceleration within the 1σ error
    # ellipse of the measurement (chi-squared 2-DOF threshold = 2.30).
    chi2_1sigma = 2.30
    f68 = mean(chi2 .<= chi2_1sigma)

    # Append consistency metrics to the text summary
    push!(stat_lines,
          @sprintf("%-20s  %10.3f  (f68=%.2f)",
                   "$(name): acc χ²_med", median(chi2), f68))

    # --- Left panel: 2D predictive scatter ---
    ax_2d = Axis(fig_acc_ppc[k, 1];
        xlabel="$(name): predicted accra [mas/yr²]",
        ylabel="$(name): predicted accdec [mas/yr²]",
        xgridvisible=false, ygridvisible=false)

    # Subsample for visual clarity; rasterize to keep file size small.
    scatter!(ax_2d, accra_pred[sample_idx], accdec_pred[sample_idx];
        color=(c, 0.3), markersize=4, rasterize=4,
        label="Posterior predictions")

    # Measured acceleration with ±1σ error bars
    errorbars!(ax_2d, [ax_meas], [ay_meas], [σ_ra];
        direction=:x, color=:black, linewidth=2)
    errorbars!(ax_2d, [ax_meas], [ay_meas], [σ_dec];
        direction=:y, color=:black, linewidth=2)
    scatter!(ax_2d, [ax_meas], [ay_meas];
        marker=:xcross, markersize=14, color=:black, strokewidth=2,
        label="Measured ± 1σ")

    k == 1 && axislegend(ax_2d; position=:rt, framevisible=false, labelsize=12)

    # --- Right panel: chi-squared residual distribution ---
    ax_chi = Axis(fig_acc_ppc[k, 2];
        xlabel="$(name): χ² (predicted vs measured, 2-DOF)",
        ylabel="Probability Density",
        xgridvisible=false, ygridvisible=false)

    hist!(ax_chi, chi2; normalization=:pdf, bins=40, color=(c, 0.7))

    # Mark the 1σ ellipse boundary.  The fraction of draws to the left of
    # this line (f₆₈) is the key consistency metric printed in the legend.
    vlines!(ax_chi, [chi2_1sigma];
        color=:black, linestyle=:dash,
        label=@sprintf("χ²=2.30 (1σ ellipse), f₆₈=%.2f", f68))
    axislegend(ax_chi; position=:rt, framevisible=false, labelsize=12)
end

save(joinpath(output_dir, "$(run_prefix)_accel_check.png"), fig_acc_ppc, px_per_unit=3)
fig_acc_ppc = nothing; GC.gc()
println("Acceleration posterior predictive check saved.")

# ── 12. IMBH position posterior + acceleration/PM overlay (two panels) ──────

println("Generating IMBH position overlay figure...")

using KernelDensity

# Smoothed KDE — 4× Silverman bandwidth to reduce contour jaggedness.
let
    n    = length(ox_samples)
    bf   = 4.2
    bw_x = 1.06 * std(ox_samples) * n^(-0.2) * bf
    bw_y = 1.06 * std(oy_samples) * n^(-0.2) * bf
    global _kde2d = kde((ox_samples, oy_samples); bandwidth = (bw_x, bw_y))
end

# σ density thresholds (1-2-3σ enclosing 68.27 / 95.45 / 99.73 % of mass).
let
    _flat = sort(vec(_kde2d.density), rev = true)
    _cum  = cumsum(_flat) ./ sum(_flat)
    global _ovl_lvls = sort([_flat[min(searchsortedfirst(_cum, f), length(_flat))]
                              for f in [0.6827, 0.9545, 0.9973]])
end

# Bounding box covering both the IMBH posterior and all star positions.
_star_ras  = [astrom_obs[n].table.ra[1]  for n in star_names]
_star_decs = [astrom_obs[n].table.dec[1] for n in star_names]
_x_lo = min(minimum(_star_ras),  quantile(ox_samples, 0.001))
_x_hi = max(maximum(_star_ras),  quantile(ox_samples, 0.999))
_y_lo = min(minimum(_star_decs), quantile(oy_samples, 0.001))
_y_hi = max(maximum(_star_decs), quantile(oy_samples, 0.999))

# Symmetric square limits centred on the origin (AvdM10 = 0,0).
_R         = max(abs(_x_lo), abs(_x_hi), abs(_y_lo), abs(_y_hi)) * 1.12
_axis_span = 2 * _R

# Adaptive arrow scales.  Acc is 50% of the original 0.20 coefficient;
# PM is 70% of the original 0.15 coefficient.
_acc_mags  = [sqrt(octo_utils.stars[n].acc_ra^2 + octo_utils.stars[n].acc_dec^2) for n in star_names]
_pm_mags   = [sqrt(pm_obs[n].table.pmra[1]^2    + pm_obs[n].table.pmdec[1]^2)    for n in star_names]
_max_acc   = maximum(_acc_mags)
_max_pm    = maximum(_pm_mags)
_scale_acc = 0.125  * _axis_span / _max_acc   # 50 % of previous 0.20, then +25 %
_scale_pm  = 0.1365 * _axis_span / _max_pm   # 70 % of previous 0.15, then +30 %

# Helper: draw KDE contours + reference markers shared by both panels.
function _draw_ovl_common!(ax)
    # Single contourf! call with four level boundaries → three filled bands.
    # Adding a ceiling just above the KDE maximum creates the third inter-level
    # band (inside 1σ) so each σ region gets its own explicit grey shade.
    # extendlow = :transparent leaves everything outside the 3σ contour unfilled.
    # The three-pass approach was incorrect: each pass's transparent "below
    # threshold" fill was overwriting the previous pass's grey fill in CairoMakie.
    contourf!(ax, _kde2d.x, _kde2d.y, _kde2d.density;
        levels   = [_ovl_lvls[1], _ovl_lvls[2], _ovl_lvls[3],
                    maximum(_kde2d.density) * 1.001],
        colormap = [RGBf(0.85, 0.85, 0.85), RGBf(0.65, 0.65, 0.65),
                    RGBf(0.45, 0.45, 0.45)],
        extendlow = :transparent)
    contour!(ax, _kde2d.x, _kde2d.y, _kde2d.density;
        levels = _ovl_lvls, color = :black, linewidth = 1.5)
    scatter!(ax, [0.0], [0.0];
        marker = '+', markersize = 44, color = :grey40)
    scatter!(ax, [median(ox_samples)], [median(oy_samples)];
        marker = :circle, color = :black, markersize = 24)
end

# Reference scale bar helper (lower-left visual corner).
# With xreversed=true: visual-left = large positive x.  The arrow is placed near
# the left edge; the text label is placed to its screen-right (smaller data-x = toward
# center), so it always extends inward and stays fully visible.
function _draw_ref_box!(ax, scale, ref_val, label_text, arrow_color,
                        shaft_w, tip_w, tip_l)
    _margin   = 0.05 * _axis_span
    _text_off = 0.030 * _axis_span   # gap between arrow shaft and label
    _ar_x     = _R - 0.18 * _axis_span   # near visual-left edge (large positive x)
    _ref_y0   = -_R + _margin
    _bar_len  = ref_val * scale
    arrows2d!(ax, [_ar_x], [_ref_y0], [0.0], [_bar_len];
        color = arrow_color, shaftwidth = shaft_w, tipwidth = tip_w, tiplength = tip_l)
    # align=(:left,:center) anchors the text's left edge at the given position and
    # extends rightward in screen space (= toward smaller x = plot interior).
    text!(ax, _ar_x - _text_off, _ref_y0 + _bar_len / 2;
        text = label_text, color = arrow_color, align = (:left, :center), fontsize = 20)
end

# Two-panel figure: left = PM, right = accelerations.
# Both panels use the same square axis limits and aspect=DataAspect() so that
# limits!() does not interfere with xreversed=true.
fig_ovl = Figure(size = (1700, 900), fontsize = 26)

_axis_kw = (
    xreversed    = true,
    aspect       = DataAspect(),
    limits       = (-_R, _R, -_R, _R),
    xgridvisible = false,
    ygridvisible = false,
)

ax_pm  = Axis(fig_ovl[1, 1];
    xlabel = Makie.rich("Δα*", Makie.subscript("IMBH"), " [mas]"),
    ylabel = Makie.rich("Δδ",  Makie.subscript("IMBH"), " [mas]"),
    _axis_kw...)

ax_acc = Axis(fig_ovl[1, 2];
    xlabel             = Makie.rich("Δα*", Makie.subscript("IMBH"), " [mas]"),
    yticklabelsvisible = false,
    ylabelvisible      = false,
    _axis_kw...)

# ── Left panel: PM vectors (no labels — legend lives in the right panel) ─
_draw_ovl_common!(ax_pm)

# 50 posterior orbit samples per star overplotted at low alpha.
let _pm_idx = sample_idx[1:2:end]   # ~50 evenly spaced samples (subset of the 100 in sample_idx)
    for name in star_names
        color = star_colors[name]
        s     = star_samples[name]
        obs_ra  = astrom_obs[name].table.ra[1]
        obs_dec = astrom_obs[name].table.dec[1]
        for idx in _pm_idx
            orb_s = Visual{KepOrbit}(;
                a=s.a[idx], e=s.e[idx], i=s.i[idx],
                ω=s.ω[idx], Ω=s.Ω[idx], tp=s.tp[idx],
                M=M_samples[idx], plx=plx_samples[idx])
            P_s = s.a[idx]^1.5 / sqrt(M_samples[idx])
            ts  = range(epoch_mjd, epoch_mjd + P_s * 365.25; length=300)
            ra_s  = [raoff(orbitsolve(orb_s, t)) + ox_samples[idx] for t in ts]
            dec_s = [decoff(orbitsolve(orb_s, t)) + oy_samples[idx] for t in ts]
            lines!(ax_pm, ra_s, dec_s; color=(color, 0.12), linewidth=0.6)
        end
    end
end

# Pass 1: star markers (before arrows so arrows render on top).
for name in star_names
    color   = star_colors[name]
    obs_ra  = astrom_obs[name].table.ra[1]
    obs_dec = astrom_obs[name].table.dec[1]
    scatter!(ax_pm, [obs_ra], [obs_dec];
        marker = '★', color = color, markersize = 28,
        strokecolor = :black, strokewidth = 0.5)
end

# Pass 2: PM arrows at front layer.
for name in star_names
    obs_ra  = astrom_obs[name].table.ra[1]
    obs_dec = astrom_obs[name].table.dec[1]
    pm_ra   = pm_obs[name].table.pmra[1]
    pm_dec  = pm_obs[name].table.pmdec[1]
    arrows2d!(ax_pm, [obs_ra], [obs_dec],
        [pm_ra * _scale_pm], [pm_dec * _scale_pm];
        color = :red, shaftwidth = 1.5, tipwidth = 8, tiplength = 8)
end

let
    # Fix reference vector at 60 km/s and derive the angular equivalent.
    _ref_pm_kms = 60.0
    _ref_pm_val = _ref_pm_kms / (ustrip(octo_utils.distance_kpc) * 4.74047)
    _pm_label   = @sprintf("%.3g mas yr⁻¹ = %d km s⁻¹", _ref_pm_val, round(Int, _ref_pm_kms))
    _draw_ref_box!(ax_pm, _scale_pm, _ref_pm_val, _pm_label, :red, 1.5, 8, 8)
end

# ── Right panel: acceleration vectors + uncertainty ellipses ─────────────
_draw_ovl_common!(ax_acc)

θs = range(0, 2π, length = 120)

# Pass 1: uncertainty ellipses + star markers — drawn before arrows so that
# the black acceleration arrows render on top.
for name in star_names
    color = star_colors[name]
    star  = octo_utils.stars[name]
    obs_ra  = astrom_obs[name].table.ra[1]
    obs_dec = astrom_obs[name].table.dec[1]
    acc_ra  = star.acc_ra;  acc_dec = star.acc_dec
    σ_ra    = star.sigma_acc_ra;  σ_dec = star.sigma_acc_dec

    tip_x = obs_ra  + acc_ra  * _scale_acc
    tip_y = obs_dec + acc_dec * _scale_acc
    for (n_sig, alpha) in [(3, 0.08), (2, 0.18), (1, 0.32)]
        ex = tip_x .+ n_sig * σ_ra  * _scale_acc .* cos.(θs)
        ey = tip_y .+ n_sig * σ_dec * _scale_acc .* sin.(θs)
        poly!(ax_acc, Point2f.(ex, ey);
            color = (color, alpha), strokecolor = (color, 0.4), strokewidth = 0.5)
    end

    scatter!(ax_acc, [obs_ra], [obs_dec];
        marker = '★', color = color, markersize = 28,
        strokecolor = :black, strokewidth = 0.5)
end

# Pass 2: acceleration arrows drawn last → frontmost layer.
for name in star_names
    star  = octo_utils.stars[name]
    obs_ra  = astrom_obs[name].table.ra[1]
    obs_dec = astrom_obs[name].table.dec[1]
    arrows2d!(ax_acc, [obs_ra], [obs_dec],
        [star.acc_ra * _scale_acc], [star.acc_dec * _scale_acc];
        color = :black, shaftwidth = 2.5, tipwidth = 12, tiplength = 12)
end

_draw_ref_box!(ax_acc, _scale_acc, round(_max_acc; sigdigits = 1),
    @sprintf("%.2g mas yr⁻²", round(_max_acc; sigdigits = 1)),
    :black, 2.5, 12, 12)

# Combined legend for both panels: framed box in the lower-right of the right panel.
# Built explicitly so that PM and acc entries can show a line + arrowhead compound swatch.
let
    _leg_entries = Any[
        LineElement(color = RGBf(0.85, 0.85, 0.85), linewidth = 8),
        LineElement(color = RGBf(0.65, 0.65, 0.65), linewidth = 8),
        LineElement(color = RGBf(0.45, 0.45, 0.45), linewidth = 8),
        MarkerElement(marker = '+',     markersize = 22, color = :grey40),
        MarkerElement(marker = :circle, markersize = 14, color = :black),
        # PM: short line + right-pointing triangle arrowhead
        [LineElement(color = :red, linewidth = 2,
                     linepoints = [Point2f(0.0, 0.5), Point2f(0.72, 0.5)]),
         MarkerElement(marker = :rtriangle, color = :red,
                       markersize = 12, points = [Point2f(0.88, 0.5)])],
        # Acc: short line + right-pointing triangle arrowhead
        [LineElement(color = :black, linewidth = 2.5,
                     linepoints = [Point2f(0.0, 0.5), Point2f(0.72, 0.5)]),
         MarkerElement(marker = :rtriangle, color = :black,
                       markersize = 14, points = [Point2f(0.88, 0.5)])],
    ]
    _leg_labels = [
        "3σ IMBH position", "2σ IMBH position", "1σ IMBH position",
        "AvdM10 centre", "IMBH (median)",
        "Proper motions", "2D accelerations",
    ]
    for name in star_names
        push!(_leg_entries, MarkerElement(marker = '★', color = star_colors[name],
                                          markersize = 14,
                                          strokecolor = :black, strokewidth = 0.5))
        push!(_leg_labels, "Star $name")
    end
    axislegend(ax_acc, _leg_entries, _leg_labels;
        position = :rb, framevisible = true, labelsize = 18)
end

# ── Secondary sky-coordinate axes (top = RA, right = Dec) ────────────────
# Convert mas offsets from AvdM10 to absolute RA/Dec in sexagesimal format.
_mas_to_ra_str(Δαstar_mas) = let
    ra_deg = octo_utils.ra_cm_deg + Δαstar_mas / (cosd(octo_utils.dec_cm_deg) * 3_600_000)
    ra_hr  = mod(ra_deg / 15, 24)
    h = floor(Int, ra_hr)
    m = floor(Int, (ra_hr - h) * 60)
    s = ((ra_hr - h) * 60 - m) * 60
    @sprintf("%dʰ%02dᵐ%05.2fˢ", h, m, s)
end
_mas_to_dec_str(Δδ_mas) = let
    dec_deg = octo_utils.dec_cm_deg + Δδ_mas / 3_600_000
    sgn = dec_deg >= 0 ? "+" : "−"
    d   = floor(Int, abs(dec_deg))
    m   = floor(Int, (abs(dec_deg) - d) * 60)
    s   = ((abs(dec_deg) - d) * 60 - m) * 60
    @sprintf("%s%d°%02d′%04.1f″", sgn, d, m, s)
end

# Choose tick positions: round to nearest 200 mas (or smaller) within ±_R.
_sky_step = let _raw = _axis_span / 4
    _p = 10^floor(log10(_raw))
    _r = _raw / _p
    (_r <= 1.5 ? _p : _r <= 3.0 ? 2*_p : _r <= 7.0 ? 5*_p : 10*_p)
end
_sky_tick_vals = let
    _start = ceil(Int, -_R / _sky_step) * _sky_step
    _stop  = floor(Int,  _R / _sky_step) * _sky_step
    Float64.(_start:_sky_step:_stop)
end
_ra_tick_labels  = [_mas_to_ra_str(v)  for v in _sky_tick_vals]
_dec_tick_labels = [_mas_to_dec_str(v) for v in _sky_tick_vals]

_sky_axis_kw = (
    limits            = (-_R, _R, -_R, _R),
    xreversed         = true,
    aspect            = DataAspect(),
    xaxisposition     = :top,
    yaxisposition     = :right,
    backgroundcolor   = :transparent,
    xgridvisible      = false,
    ygridvisible      = false,
    bottomspinevisible = false,
    leftspinevisible   = false,
    xticks            = (_sky_tick_vals, _ra_tick_labels),
    yticks            = (_sky_tick_vals, _dec_tick_labels),
    xticklabelsize    = 18,
    yticklabelsize    = 18,
    xlabelsize        = 20,
    ylabelsize        = 20,
)

# PM panel: top RA axis only (Dec labels would duplicate the acc panel).
Axis(fig_ovl[1, 1];
    xlabel             = "α (J2000)",
    yticklabelsvisible = false,
    ylabelvisible      = false,
    rightspinevisible  = false,
    _sky_axis_kw...)

# Acc panel: top RA + right Dec axes.
Axis(fig_ovl[1, 2];
    xlabel = "α (J2000)",
    ylabel = "δ (J2000)",
    _sky_axis_kw...)

save(joinpath(output_dir, "$(run_prefix)_imbh_position.png"), fig_ovl, px_per_unit = 3)
fig_ovl = nothing; GC.gc()
println("IMBH position overlay figure saved.")

# ── 13. 3D orbit animation (360° pan, IMBH-centric) ───────────────────────

println("Generating 3D orbit animation...")

using PlanetOrbits: posx, posy, posz

const AU_PER_PC = 206265.0

# Pre-compute 3D orbit trajectories for each posterior sample (reuses
# sample_idx from the orbit panels).  All positions are relative to the
# IMBH and converted from AU to parsec.
orbit_3d_samples = Dict{String, Vector{NamedTuple}}()
star_pos_3d      = Dict{String, Vector{NamedTuple}}()

for name in star_names
    s = star_samples[name]
    orbits_list = NamedTuple[]
    pos_list    = NamedTuple[]
    for idx in sample_idx
        orb = Visual{KepOrbit}(;
            a = s.a[idx], e = s.e[idx], i = s.i[idx],
            ω = s.ω[idx], Ω = s.Ω[idx], tp = s.tp[idx],
            M = M_samples[idx], plx = plx_samples[idx])
        P_yr = s.a[idx]^1.5 / sqrt(M_samples[idx])
        ts   = range(epoch_mjd, epoch_mjd + P_yr * 365.25; length=300)
        sols = [orbitsolve(orb, t) for t in ts]
        push!(orbits_list, (
            x = [posx(sl) for sl in sols] ./ AU_PER_PC,
            y = [posy(sl) for sl in sols] ./ AU_PER_PC,
            z = [posz(sl) for sl in sols] ./ AU_PER_PC,
        ))
        sol_now = orbitsolve(orb, epoch_mjd)
        push!(pos_list, (
            x = posx(sol_now) / AU_PER_PC,
            y = posy(sol_now) / AU_PER_PC,
            z = posz(sol_now) / AU_PER_PC,
        ))
    end
    orbit_3d_samples[name] = orbits_list
    star_pos_3d[name]      = pos_list
end

# Axis limit matched to the combined 2D orbit panel (half_all in mas → pc)
# raoff [mas] = posx [AU] × plx [mas], so half_all [mas] / (plx [mas] × AU_PER_PC) = half [pc]
lim = half_all / (median(plx_samples) * AU_PER_PC) * 0.7

# Animation view-angle parameters (also used for the initial Axis3 view).
azim_start = -3 * π / 4    # starting azimuth (three-quarter view)
elev_max   = 50 * π / 180  # start/end elevation
elev_min   = 10 * π / 180  # low angle, reveals LOS depth

# Build figure with Axis3.
# viewmode = :fit keeps ticks/labels inside the viewport as the camera rotates
# (the default :fitzoom lets labels drift outside the frame).
# dark=true: black background + white labels + white IMBH border (for animation and dark still).
# dark=false: default theme (for the light-theme still frame).
function _build_3d_fig(dark::Bool)
    fig = Figure(size = (800, 800), fontsize = 16, figure_padding = 30)
    ax  = Axis3(fig[1, 1];
        xlabel = "x [pc]", ylabel = "y [pc]", zlabel = "z (LOS) [pc]",
        limits = (-lim, lim, -lim, lim, -lim, lim),
        aspect = :data,
        viewmode = :fit,
        azimuth   = azim_start,
        elevation = elev_max,
    )

    sample_masses_loc = M_samples[sample_idx]
    mass_ref_loc      = median(sample_masses_loc)
    base_size_loc     = lim * 0.012

    for (k, name) in enumerate(star_names)
        color = star_colors[name]
        for o in orbit_3d_samples[name]
            lines!(ax, o.x, o.y, o.z; color = (color, 0.3), linewidth = 0.5)
        end
        px = [p.x for p in star_pos_3d[name]]
        py = [p.y for p in star_pos_3d[name]]
        pz = [p.z for p in star_pos_3d[name]]
        marker_sizes = base_size_loc .* (sample_masses_loc ./ mass_ref_loc)
        meshscatter!(ax, px, py, pz; markersize = marker_sizes, color = color)
        lines!(ax, [NaN], [NaN], [NaN]; color = color, linewidth = 2, label = "Star $name")
    end

    # IMBH at origin.  On a dark background the black sphere needs a white border
    # for readability: draw a slightly larger white sphere first, then the black one.
    if dark
        meshscatter!(ax, [0.0], [0.0], [0.0]; markersize = lim * 0.020, color = :white)
    end
    meshscatter!(ax, [0.0], [0.0], [0.0]; markersize = lim * 0.015, color = :black, label = "IMBH")

    axislegend(ax; position = :rt, framevisible = false)
    return fig, ax
end

n_frames   = 240
framerate  = 15
anim_path  = joinpath(output_dir, "$(run_prefix)_orbits_3d.mp4")

# Still-frame view: ≈7.5 s into the animation (frame 112 at 15 fps).
_still_frame  = 112
_t_still      = _still_frame / n_frames
_azim_still   = azim_start + 2π * _t_still
_elev_still   = (elev_max + elev_min) / 2 +
                (elev_max - elev_min) / 2 * cos(2π * _t_still)

# Dark-theme animation + dark-theme still frame.
with_theme(theme_dark()) do
    fig3d, ax3 = _build_3d_fig(true)

    record(fig3d, anim_path, 0:(n_frames - 1); framerate) do frame
        t = frame / n_frames
        ax3.azimuth[]   = azim_start + 2π * t
        ax3.elevation[] = (elev_max + elev_min) / 2 +
                          (elev_max - elev_min) / 2 * cos(2π * t)
    end
    println("3D orbit animation saved to: $anim_path")

    ax3.azimuth[]   = _azim_still
    ax3.elevation[] = _elev_still
    still_dark_path = joinpath(output_dir, "$(run_prefix)_orbits_3d_still_dark.png")
    save(still_dark_path, fig3d, px_per_unit = 3)
    println("3D still frame (dark) saved to: $still_dark_path")
end

# Light-theme animation + still frame (default Makie theme).
let
    anim_light_path = joinpath(output_dir, "$(run_prefix)_orbits_3d_light.mp4")
    fig3d_lt, ax3_lt = _build_3d_fig(false)
    record(fig3d_lt, anim_light_path, 0:(n_frames - 1); framerate) do frame
        t = frame / n_frames
        ax3_lt.azimuth[]   = azim_start + 2π * t
        ax3_lt.elevation[] = (elev_max + elev_min) / 2 +
                              (elev_max - elev_min) / 2 * cos(2π * t)
    end
    println("3D orbit animation (light) saved to: $anim_light_path")
    ax3_lt.azimuth[]   = _azim_still
    ax3_lt.elevation[] = _elev_still
    still_light_path = joinpath(output_dir, "$(run_prefix)_orbits_3d_still.png")
    save(still_light_path, fig3d_lt, px_per_unit = 3)
    println("3D still frame (light) saved to: $still_light_path")
end

# ── 14. Observed vs. model residuals (position, PM, acceleration) ────────────

println("Generating residuals figure...")

# For each posterior draw in sample_idx, evaluate all 6 kinematic model predictions
# at epoch_mjd:
#   pos (Δα*, Δδ)       [mas]      — orbital raoff/decoff + IMBH offset
#   PM  (μ_α*, μ_δ)     [mas/yr]   — central finite difference (half-step _dt_fd = 0.5 yr)
#   acc (a_α*, a_δ)     [mas/yr²]  — PlanetOrbits accra/accdec
# Pull = (observed − median_prediction) / σ_obs.
# Error bar on each point = std(prediction) / σ_obs (posterior predictive spread).

_dt_fd = 0.5   # finite-difference half-interval [yr]

_n_obs_res = 6
_pull      = zeros(_n_obs_res, length(star_names))
_pstd      = zeros(_n_obs_res, length(star_names))

for (k, name) in enumerate(star_names)
    s   = star_samples[name]
    n_s = length(sample_idx)
    preds = zeros(n_s, _n_obs_res)

    for (j, idx) in enumerate(sample_idx)
        orb   = Visual{KepOrbit}(;
            a=s.a[idx], e=s.e[idx], i=s.i[idx],
            ω=s.ω[idx], Ω=s.Ω[idx], tp=s.tp[idx],
            M=M_samples[idx], plx=plx_samples[idx])
        sol0  = orbitsolve(orb, epoch_mjd)
        sol_p = orbitsolve(orb, epoch_mjd + _dt_fd * 365.25)
        sol_m = orbitsolve(orb, epoch_mjd - _dt_fd * 365.25)

        preds[j, 1] = raoff(sol0)  + ox_samples[idx]
        preds[j, 2] = decoff(sol0) + oy_samples[idx]
        preds[j, 3] = (raoff(sol_p)  - raoff(sol_m))  / (2 * _dt_fd)
        preds[j, 4] = (decoff(sol_p) - decoff(sol_m)) / (2 * _dt_fd)
        preds[j, 5] = accra(sol0)
        preds[j, 6] = accdec(sol0)
    end

    obs_vals = [
        astrom_obs[name].table.ra[1],   astrom_obs[name].table.dec[1],
        pm_obs[name].table.pmra[1],     pm_obs[name].table.pmdec[1],
        acc_obs[name].table.accra[1],   acc_obs[name].table.accdec[1],
    ]
    obs_sigs = [
        astrom_obs[name].table.σ_ra[1],   astrom_obs[name].table.σ_dec[1],
        pm_obs[name].table.σ_pmra[1],     pm_obs[name].table.σ_pmdec[1],
        acc_obs[name].table.σ_accra[1],   acc_obs[name].table.σ_accdec[1],
    ]

    for i in 1:_n_obs_res
        _pull[i, k] = (obs_vals[i] - median(preds[:, i])) / obs_sigs[i]
        _pstd[i, k] = std(preds[:, i]) / obs_sigs[i]
    end
end

_res_ylabels = [
    Makie.rich("Δα* residual [σ]"),
    Makie.rich("Δδ residual [σ]"),
    Makie.rich("μ", Makie.subscript("α*"), " residual [σ]"),
    Makie.rich("μ", Makie.subscript("δ"),  " residual [σ]"),
    Makie.rich("a", Makie.subscript("α*"), " residual [σ]"),
    Makie.rich("a", Makie.subscript("δ"),  " residual [σ]"),
]

fig_res = Figure(size = (900, 700), fontsize = 18)
_xs_res = Float64.(1:length(star_names))

for i in 1:_n_obs_res
    row = ceil(Int, i / 2)
    col = mod1(i, 2)
    ax  = Axis(fig_res[row, col];
        ylabel       = _res_ylabels[i],
        xticks       = (_xs_res, star_names),
        xgridvisible = false,
        ygridvisible = false,
    )
    hspan!(ax, -1.0, 1.0; color = (:grey, 0.12))
    hlines!(ax, [0.0]; color = :grey40, linestyle = :dash, linewidth = 1)

    for (k, name) in enumerate(star_names)
        c = star_colors[name]
        errorbars!(ax, [_xs_res[k]], [_pull[i, k]], [_pstd[i, k]];
            color = (c, 0.5), linewidth = 2, whiskerwidth = 6)
        scatter!(ax, [_xs_res[k]], [_pull[i, k]];
            color = c, markersize = 10, strokecolor = :black, strokewidth = 0.5)
    end

    xlims!(ax, 0.5, length(star_names) + 0.5)
    row < 3 && (ax.xticklabelsvisible[] = false)
end

save(joinpath(output_dir, "$(run_prefix)_residuals.png"), fig_res, px_per_unit = 3)
fig_res = nothing; GC.gc()
println("Residuals figure saved.")

# ── 15. Write collected posterior stats (after all diagnostic sections) ──────

open(stats_path, "w") do io
    println(io, "Posterior summaries (median, 68% CI)")
    println(io, "Chain: $chain_path")
    println(io, "Epoch: $epoch_mjd MJD ($epoch_year yr)")
    println(io, "Stars: $(join(star_names, ", "))")
    println(io)
    for line in stat_lines
        println(io, line)
    end
end
println("Posterior stats saved to: $stats_path")

println("\nDone. All plots written to: $output_dir")
