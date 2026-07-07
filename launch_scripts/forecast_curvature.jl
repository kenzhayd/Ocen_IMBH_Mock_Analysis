"""
Forecast orbital curvature for future observations.

Uses chain from the primary real-data fit to compute expected 
growth of curvature residuals over time. 
Suggests a timeline for when orbital curvature will 
be significant enough to constrain accelerations.
"""

ENV["OCTOFITTERPY_AUTOLOAD_EXTENSIONS"] = "yes"

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
using DataFrames
using Printf
using Random
using CairoMakie
using PairPlots
using LaTeXStrings


push!(LOAD_PATH, @__DIR__)
include("octo_utils.jl")

# ===================== CONFIGURATION =====================

# All outputs are saved to output_dir
output_dir = raw"C:\Users\macke\Clusters\Ocen_IMBH_Mock_Analysis\mock_results\forcast_curvature\updated_5000_samples_2"
mkpath(output_dir) 

# Number of posterior orbit samples used to estimate the
# distribution of future curvature residuals
N_samp = 5000

# Number of posterior orbits drawn on the orbit panel. This matches the
# plot_chain.jl sky-plane orbit panels.
N_orbit_plot = 50

# Assumed astrometric uncertainty of position measurements (mas)
sigma_astrometry = 0.5

# Reference epoch (Same as chain)
epoch_year = 2010.0
epoch_mjd  = Octofitter.years2mjd(epoch_year)

# Same star colour convention as plot_chain.jl's orbit panels.
star_names = ["A", "C", "D", "E", "F"]
w = Makie.wong_colors()
star_colors = Dict{String, Any}(
    name => w[mod1(k + 1, length(w))]
    for (k, name) in enumerate(star_names)
)

# Years to forcast orbital curvature
start_year = 2023.0
end_year   = 2035.0

forecast_years = start_year:end_year

# Time elapsed since the reference epoch 
years_since_obs = forecast_years .- epoch_year

# ===================== CHAIN =====================

# Posterior chain from real-data Octofitter run
chain_path = raw"C:\Users\macke\Clusters\run_108836842\starsACDEF_192c_18r_cont_10836842_chain.fits"
chain = Octofitter.loadchain(chain_path)

# Load chain samples
logpost = vec(Array(chain[:logpost]))
idx_all = eachindex(logpost)


# ===================== FUNCTIONS =====================


"""
Extract a single posterior sample from the chain and construct the
corresponding Keplerian orbit.

Given a chain index and star name, gets the orbital parameters,
computes the semimajor axis and time of periastron passage, and returns both
the resulting orbit object and its parameters.

Returns
-------
orbit
    A Visual{KepOrbit} object.
sample
    Named tuple containing the sampled orbital and astrometric parameters.
"""
function extract_orbit_sample(chain, idx, star, epoch_mjd)

    sample = (
        P = chain[Symbol("$(star)_P")][idx],
        e = chain[Symbol("$(star)_e")][idx],
        i = chain[Symbol("$(star)_i")][idx],
        ω = chain[Symbol("$(star)_ω")][idx],
        Ω = chain[Symbol("$(star)_Ω")][idx],
        θ = chain[Symbol("$(star)_θ")][idx],

        M = chain[:M][idx],
        plx = chain[:plx][idx],

        offsetx = chain[:offsetx][idx],
        offsety = chain[:offsety][idx]
    )

    a = cbrt(sample.M * sample.P^2)

    tp = θ_at_epoch_to_tperi(
        sample.θ,
        epoch_mjd;
        a=a,
        e=sample.e,
        i=sample.i,
        ω=sample.ω,
        Ω=sample.Ω,
        M=sample.M
    )

    orbit = Visual{KepOrbit}(;
        a,
        e=sample.e,
        i=sample.i,
        ω=sample.ω,
        Ω=sample.Ω,
        tp,
        M=sample.M,
        plx=sample.plx
    )

    return orbit, sample
end


"""
Build orbit trajectories from randomly selected posterior samples.

Draws N samples from the posterior chain, propagates each orbit over one
orbital period, and returns the corresponding RA and Dec offsets.

Returns
-------
ra_samples
    Vector of RA offsets in mas.
dec_samples
    Vector of Dec offsets in mas.
"""
function sample_orbit_trajectory(chain, idx_all, star, epoch_mjd, N)

    ra_samples  = Vector{Vector{Float64}}(undef, N)
    dec_samples = Vector{Vector{Float64}}(undef, N)

    for k in 1:N
        idx = rand(idx_all)

        orbit, sample = extract_orbit_sample(chain, idx, star, epoch_mjd)

        ts = range(
            epoch_mjd,
            epoch_mjd + sample.P * 365.25,
            length = 400
        )

        sols = orbitsolve.(orbit, ts)

        ra_samples[k]  = raoff.(sols) .+ sample.offsetx
        dec_samples[k] = decoff.(sols) .+ sample.offsety
    end

    return ra_samples, dec_samples
end

"""
Compute the distribution of curvature residuals relative to a linear model.

For each forecast epoch, draws N_samp posterior orbit samples, 
builds the corresponding orbits, 
and compares them to a constant proper-motion trajectory.

Returns
-------
pos_residual
    Total positional residuals (mas).
ra_residuals
    RA residuals (mas).
dec_residuals
    Dec residuals (mas).
"""

function curvature_residuals_distribution(chain, idx_all, starname, years_since_obs, epoch_mjd; N_samp)

    pos_residual = zeros(length(years_since_obs), N_samp)
    ra_residuals = zeros(length(years_since_obs), N_samp)
    dec_residuals = zeros(length(years_since_obs), N_samp)

    star = octo_utils.stars[starname]

    # Reference epoch positions (using observed absolute positions) relative from cluster center
    ra0 = (star.ra - octo_utils.ra_cm_deg) * 3600 * 1000 * cosd(octo_utils.dec_cm_deg)
    dec0 = (star.dec - octo_utils.dec_cm_deg) * 3600 * 1000

    # Linear trajectory using observed pm and positions
    ra_lin  = ra0 .+ star.pm_ra .* years_since_obs
    dec_lin = dec0 .+ star.pm_dec .* years_since_obs

    times = epoch_mjd .+ years_since_obs .* 365.25

    for k in 1:N_samp
        idx = rand(idx_all)
        orbit, sample = extract_orbit_sample(chain, idx, starname, epoch_mjd)

        sols = orbitsolve.(orbit, times)

        # Switch to cluster-centered coordinate system from IMBH-centered
        ra_orbit  = raoff.(sols) .+ sample.offsetx
        dec_orbit = decoff.(sols) .+ sample.offsety

        # Compute residuals (Octofitter Model - Linear Model (from observed data))
        ra_residuals[:, k]  = ra_orbit .- ra_lin
        dec_residuals[:, k] = dec_orbit .- dec_lin

        pos_residual[:, k] = hypot.(ra_residuals[:, k], dec_residuals[:, k])
    end

    return pos_residual, ra_residuals, dec_residuals
end

# ===================== CURVATURE FORCAST ANALYSIS =====================

# Open a file to store curbature forcast stats for all stars
stats_file_path = joinpath(output_dir, "curvature_forecast_stats.txt")
open(stats_file_path, "w") do stats_io 
    write(stats_io, """
    Chain: $chain_path
    Epoch: $epoch_year yr ($epoch_mjd MJD)
    N_samp: $N_samp
    N_orbit_plot: $N_orbit_plot
    Sigma_astrometry: $sigma_astrometry mas
    Forecast Years: $start_year to $end_year
    """)
    write(stats_io, "\n" * "="^50 * "\n") 

frac_3sigma_by_star = Dict{String, Vector{Float64}}()

# Main loop 
for star in star_names

    # ================= OCTOFITTER MODEL POSITIONS =================
    offsetx = median(chain[:offsetx])
    offsety = median(chain[:offsety])

    star_color = star_colors[star]
    ra_samples, dec_samples = sample_orbit_trajectory(chain, idx_all, star, epoch_mjd, N_orbit_plot)
    
    # ================= OBSERVED POSITIONS =================
    real = octo_utils.stars[star]

    obs_ra  = (real.ra  - octo_utils.ra_cm_deg)  * 3600 * 1000 * cosd(octo_utils.dec_cm_deg)
    obs_dec = (real.dec - octo_utils.dec_cm_deg) * 3600 * 1000

    # ================= RESIDUALS =================
   
    pos_residuals, ra_residuals, dec_residuals = curvature_residuals_distribution(
    chain,
    idx_all,
    star,
    years_since_obs,
    epoch_mjd;
    N_samp
    )
    
    # Median of residuals for each year forcasted
    med_residual = mapslices(median, pos_residuals; dims=2)[:]

    # 1-Sigma spread (16th and 84th percentiles) for each year forcasted
    residual_lo  = mapslices(x -> quantile(x, 0.16), pos_residuals; dims=2)[:]
    residual_hi  = mapslices(x -> quantile(x, 0.84), pos_residuals; dims=2)[:]

    # 3-Sigma spread (0.15th and 99.85th percentiles) for each year forcasted
    three_sigma_lo = mapslices(x -> quantile(x, 0.0015), pos_residuals; dims=2)[:]
    three_sigma_hi = mapslices(x -> quantile(x, 0.9985), pos_residuals; dims=2)[:]

    # Fraction of residuals above 1-sigma of astrometric uncertainty for each year
    frac_1sigma = mapslices(x -> sum(x .> sigma_astrometry) / length(x), pos_residuals; dims=2)[:]

    # Fraction of residuals above 3-sigma of astrometric uncertainty for each year
    frac_3sigma = mapslices(x -> sum(x .> 3 * sigma_astrometry) / length(x), pos_residuals; dims=2)[:]
    frac_3sigma_by_star[star] = frac_3sigma
    
    # Linear Trajectory (2002 to end_year) 
    baseline_years = 2002.0:0.1:end_year
    years_since_baseline_epoch = baseline_years .- epoch_year  

    # Full linear trajectory from baseline epoch to end_year
    ra_full = obs_ra .+ real.pm_ra .* years_since_baseline_epoch 
    dec_full = obs_dec .+ real.pm_dec .* years_since_baseline_epoch 
   

    # Save residual stats 
        io_buffer = IOBuffer()
        println(io_buffer, "\n--- Star $star ---")
        println(io_buffer, "Curvature residual growth (|Δr|, mas):")
        println(io_buffer, "  Year    Median  ±1σ (68% CI)        ±3σ (99.7% CI)        N(|Δr|>σ)  N(|Δr|>3σ)"),    
        println(io_buffer, "          [mas]   [+/- mas]           [+/- mas]             /N_tot]    /N_tot")

        for (yr, m, lo1, hi1, lo3, hi3, f1, f3) in zip(
            forecast_years,
            med_residual,
            residual_lo,
            residual_hi,
            three_sigma_lo,
            three_sigma_hi,
            frac_1sigma,
            frac_3sigma
        )
            @printf(io_buffer,
                "  %5d :%6.3f   (+%.3f / -%.3f)   (+%.3f / -%.3f)     %.3f      %.3f\n",
                yr,
                m,
                hi1 - m, m - lo1,
                hi3 - m, m - lo3,
                f1, f3
            )
        end

        stats_string = String(take!(io_buffer))

        print(stats_string)
        write(stats_io, stats_string)


    # ================= FIGURE 1: ORBITS AND RESIDUAL GROWTH =================
  
    fig = Figure(size = (900, 800), 
                 layout = GridLayout(2, 2))
    

    ax1 = Axis(fig[1, 1];
        xlabel="Δα* [mas]",
        ylabel="Δδ [mas]",
        xreversed=true,
        titlealign =:left,
        ylabelsize = 26,
        xlabelsize = 26,
        xticklabelsize = 24,
        yticklabelsize = 24
        
    )

    ax2 = Axis(fig[2, 1:2];
    xlabel="Year",
    ylabel="Residual [mas]",
    titlealign =:left,
    xticks = [2025, 2030, 2035],
    xminorticks = 2023:1:2035,
    xminorticksvisible=true,
    xminorticksize=6,
    ylabelsize = 26,
    xlabelsize = 26,
    xticklabelsize = 24,
    yticklabelsize = 24
    )

    Colorbar(fig[1, 2];
         colormap = :greys,
         limits = (minimum(baseline_years), maximum(baseline_years)),
         width = 25,
         height = Auto(),
         ticks = [2005, 2015, 2025, 2035],
         minorticks = 2002:1:2035,
         minorticksvisible=false,
         minorticksize=6,
         ticklabelsize = 24,
         labelrotation = π/2,   
         valign = :center,
         halign = :left)

    # ================= ORBIT PLOTS =================
    
    # Plot orbit trajectories
    for k in eachindex(ra_samples)
        ra = ra_samples[k]
        dec = dec_samples[k]

        lines!(
            ax1,
            ra,
            dec;
            color = (star_color, 0.5),
            linewidth = 0.5,
            alpha = 0.9
        )
    end
    

    # Plot timeline colourbar 
    lines!(ax1,
        ra_full,
        dec_full;
        color = baseline_years,  
        colormap = :greys,       
        linewidth = 8,      
    )   
    
    # Overlay observed baseline
    obs_mask = baseline_years .<= 2023

    lines!(ax1,
        ra_full[obs_mask],
        dec_full[obs_mask];
        color = :white,
        linewidth = 2,
        linestyle = :dash,
    )

    # Calculate limits 
    ra_min, ra_max = extrema(ra_full)
    dec_min, dec_max = extrema(dec_full)

    # Apply calculated limits to zoom the plot
    xlims!(ax1, ra_min, ra_max)
    ylims!(ax1, dec_min, dec_max)

    # Plot observed star position at the reference epoch
    scatter!(ax1, [obs_ra], [obs_dec];
        color= w[2], markersize=18, marker=:star5, 
        strokecolor=:black, strokewidth=0.5, label= "Star $star"
    )
    
    # Star labels
    if star in ["E", "D", "F"]
        text!(ax1, "Star $star";
            position=(0.97, 0.95),
            align=(:right, :top),
            space=:relative,
            fontsize=24
        )
    else 
        text!(ax1, "Star $star"; position=(0.02, 0.95), align=(:left, :top),
        space=:relative, fontsize=24)
    end

    # ================= RESIDUAL DISTRIBUTION PLOT =================
    
    # Plot residuals for total position (Octofitter trajectory - linear trajectory)
    lines!(ax2, forecast_years, med_residual;
        color = star_color,
        linewidth = 3,
        label = "Median |Δr|"
    )

    # 1-Sigma spread
    band!(ax2, forecast_years, residual_lo, residual_hi;
    color = star_color, alpha=0.2, label="±1σ |Δr|")

    
    # Astrometric 3-sigma uncertainty
    hlines!(ax2, 3*sigma_astrometry;
        color=:black, linestyle=:dot, label="3σ detection",
        linewidth=2
    )

    ylims!(ax2, 0, nothing)
    xlims!(ax2, start_year, end_year)
    axislegend(ax2; position=:lt, framevisible=true, backgroundcolor=:white, labelsize=24)

    display(fig)

    # Save plot to output_dir
    plot_filename = joinpath(output_dir, "curvature_forecast_star$(star).png")
        save(plot_filename, fig)
        println("Saved plot for Star $star to: $plot_filename")

end

# ================= FIGURE 2: FRACTION ABOVE 3-SIGMA  =================

fig2 = Figure(size = (650, 550),
            figure_padding = (10, 30, 10, 20))  # left, right, bottom, top

ax_frac = Axis(fig2[1, 1];
    xlabel="Year",
    ylabel="3σ Detection Probability",
    xticks = [2025,2030, 2035],
    xminorticks = 2023:1:2035,
    xminorticksvisible=true,
    xminorticksize=6,
    yticks = [0.0, 0.25, 0.5, 0.75, 1.0],
    ylabelsize = 26,
    xlabelsize = 26,
    xticklabelsize = 24,
    yticklabelsize = 24
)

for star in star_names
    lines!(ax_frac, forecast_years, frac_3sigma_by_star[star];
        color = star_colors[star],
        linewidth = 3,
        label = "Star $star"
    )
end

ylims!(ax_frac, 0, 1)
xlims!(ax_frac, start_year, end_year)
axislegend(ax_frac; position=:lt, framevisible=true, backgroundcolor=:white, nbanks=3, fontsize =20)

display(fig2)

plot_filename_frac = joinpath(output_dir, "residual_fractions_3sigma_all_stars.png")
save(plot_filename_frac, fig2)
println("Saved combined 3-sigma fractions plot to: $plot_filename_frac")

end 
println("\nAll saved to directory: $output_dir")




