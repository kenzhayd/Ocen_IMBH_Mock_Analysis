"""
Batches of posterior plots.

Generates posterior histogram panels with input/recovered median lines.
Reuses plotting code from plot_chain.jl (Section 11).
"""

using TOML
using CairoMakie
using Statistics
using Printf
using Octofitter
using Distributions

push!(LOAD_PATH, @__DIR__)
include(joinpath(@__DIR__,"parse_config.jl"))

# ===================== CONFIGURATION =====================

# Usage: julia batch_posterior_plots.jl <directory> OR set run_directory variable.

# Direcotry with runs to plot 
run_directory = raw"C:\Users\macke\Clusters\Ocen_IMBH_Mock_Analysis\mock_results_lite\rand2_192c_16r_June20"
stars = ["A", "C", "D", "E", "F"] 

# Output subdirectory
output_dir = joinpath(run_directory, "batch_plots")
mkpath(output_dir)

# ===================== FUNCTIONS =====================

"""
Extract TOML block from summary markdown file. 
Reused from plot_chain.jl.
"""
function extract_toml_from_summary(summary_text::String)
    m = match(r"`toml\r?\n(.*?)`"s, summary_text)
    m === nothing && error("Could not find a TOML block in summary file")
    return TOML.parse(m[1])
end

"""
Extract input parameters from mock config.
"""
function get_input_parameters(cfg::Dict,
                              epoch_mjd::Float64,
                              stars)
    mock = get(cfg, "mock", Dict())
    params = Dict{String, Float64}()
    
    # System params
    params["M"] = get(mock, "M_IMBH", NaN)
    params["plx"] = get(mock, "plx", NaN)
    params["offsetx"] = get(mock, "offsetx", NaN)
    params["offsety"] = get(mock, "offsety", NaN)
    
    # Star params
    stars_cfg = get(mock, "stars", Dict())
    for name in stars
        scfg = get(stars_cfg, name, Dict())
        oe = get(scfg, "orbital_elements", Float64[])
        
        P, e, i, ω, Ω, θ = oe[1:6]
        M = params["M"]
        a = cbrt(M * P^2)
        
        params["$(name)_a"] = a
        params["$(name)_e"] = e
        params["$(name)_i"] = i      
        params["$(name)_ω"] = ω     
        params["$(name)_Ω"] = Ω      
    end
    return params
end

"""
Adapted param_panel! function from plot_chain.jl (Section 11).
Added vertical lines for input truth and fitted median.
"""
function param_panel!(layout, row, col, color, samples, xlabel, output_dir, run_prefix, fig_post;
                      input_val=nothing, median_val=nothing,
                      show_legend=false, xlims=nothing, bins=30, xticks=Makie.automatic)
    ax = Axis(layout[row, col]; xlabel=xlabel, ylabel="Probability Density",
              xgridvisible=false, ygridvisible=false, xticks=xticks, width=230, height=220)
    # Plot histogram
    hist!(ax, samples; normalization=:pdf, bins=bins, color=(color, 0.7))
    # Median recovered
    if median_val !== nothing
        vlines!(ax, [median_val]; color=Makie.wong_colors()[2], linestyle=:solid, label="Recovered Median")
    end
    # Input parameter
    if input_val !== nothing
        vlines!(ax, [input_val]; color=Makie.wong_colors()[3], linestyle=:solid, label="Input")
    end

    show_legend && axislegend(ax; position=:lt, framevisible=false)
    xlims !== nothing && Makie.xlims!(ax, xlims...)
end


"""
Generate posterior plot for a single run.
Adapted code from plot_chain.jl (Section 11).
"""
function generate_posterior_plot(run_dir, run_prefix, input_vals, stars, output_dir)
    
    chain_path = joinpath(run_dir, "$(run_prefix)_chain.fits")
    isfile(chain_path) || error("Chain file not found: $chain_path")
    
    println("  Loading chain: $chain_path")
    chain = Octofitter.loadchain(chain_path)

    # Extract samples 
    M_samples = vec(chain[:M])
    plx_samples = vec(chain[:plx])
    ox_samples = vec(chain[:offsetx])
    oy_samples = vec(chain[:offsety])

    star_samples = Dict{String, NamedTuple}()
    for name in stars
    star_samples[name] = (
        a = vec(chain[Symbol("$(name)_a")]),
        e = vec(chain[Symbol("$(name)_e")]),
        i = vec(chain[Symbol("$(name)_i")]),
        ω = vec(chain[Symbol("$(name)_ω")]),
        Ω = vec(chain[Symbol("$(name)_Ω")]),
    )
    end

    # Compute medians 
    medians = Dict{String, Float64}(
        "M" => median(M_samples),
        "plx" => median(plx_samples),
        "offsetx" => median(ox_samples),
        "offsety" => median(oy_samples),
    )
    for name in stars
        s = star_samples[name]
        medians["$(name)_a"]  = median(s.a)
        medians["$(name)_e"]  = median(s.e)
        medians["$(name)_i"]  = median(s.i)
        medians["$(name)_ω"]  = median(s.ω)
        medians["$(name)_Ω"]  = median(s.Ω)
    end
    
    # Colors 
    sys_color = Makie.wong_colors()[1]
    star_colors = Dict(name => Makie.wong_colors()[mod1(k+1, length(Makie.wong_colors()))] 
                       for (k, name) in enumerate(stars))

    # Build figure 
    n_stars = length(stars)
    fig_post = Figure(size=(1800, (1 + n_stars) * 300), fontsize=18)

    iv(key) = get(input_vals, key, nothing)
    mv(key) = get(medians, key, nothing)
    
    # System row 
    param_panel!(fig_post, 1, 1, sys_color, M_samples ./ 1e4,
                 Makie.rich("M", Makie.subscript("IMBH"), " [10 M", Makie.subscript("☉"), "]"),
                 output_dir, run_prefix, fig_post;
                 input_val=iv("M")/1e4, median_val=mv("M")/1e4)
                 
    param_panel!(fig_post, 1, 2, sys_color, plx_samples, "plx [mas]",
                 output_dir, run_prefix, fig_post;
                 input_val=iv("plx"), median_val=mv("plx"))
                 
    param_panel!(fig_post, 1, 3, sys_color, ox_samples,
                 Makie.rich("Δα*", Makie.subscript("IMBH"), " [mas]"),
                 output_dir, run_prefix, fig_post;
                 input_val=iv("offsetx"), median_val=mv("offsetx"))
                 
    param_panel!(fig_post, 1, 4, sys_color, oy_samples,
                 Makie.rich("Δδ", Makie.subscript("IMBH"), " [mas]"),
                 output_dir, run_prefix, fig_post;
                 input_val=iv("offsety"), median_val=mv("offsety"))

    lines = [
    LineElement(color=Makie.wong_colors()[2]),
    LineElement(color=Makie.wong_colors()[3])
    ]
    labels = ["Recovered Median", "Input Value"]

    Legend(fig_post[1, 5], lines, labels; framevisible=true, labelsize=16)

    for (k, name) in enumerate(stars)
        row = k + 1
        c = star_colors[name]
        s = star_samples[name]
        
        # Convert angles
        i_deg   = rad2deg.(s.i)
        ω_deg   = rad2deg.(s.ω)
        Ω_deg   = rad2deg.(s.Ω)
        iv_i    = iv("$(name)_i")  !== nothing ? rad2deg(iv("$(name)_i"))  : nothing
        iv_ω    = iv("$(name)_ω")  !== nothing ? rad2deg(iv("$(name)_ω"))  : nothing
        iv_Ω    = iv("$(name)_Ω")  !== nothing ? rad2deg(iv("$(name)_Ω"))  : nothing
        mv_i    = rad2deg(mv("$(name)_i"))
        mv_ω    = rad2deg(mv("$(name)_ω"))
        mv_Ω    = rad2deg(mv("$(name)_Ω"))
       
        param_panel!(fig_post, row, 1, c, s.a, "$(name): a [AU]",
                     output_dir, run_prefix, fig_post;
                     input_val=iv("$(name)_a"), median_val=mv("$(name)_a"),
                     xlims=(0, 20_000), bins=range(0, 20_000, length=51),
                     xticks=[0, 10_000, 20_000])
        param_panel!(fig_post, row, 2, c, s.e, "$(name): e",
                     output_dir, run_prefix, fig_post;
                     input_val=iv("$(name)_e"), median_val=mv("$(name)_e"))
        param_panel!(fig_post, row, 3, c, i_deg, "$(name): i [°]",
                     output_dir, run_prefix, fig_post;
                     input_val=iv_i, median_val=mv_i)
        param_panel!(fig_post, row, 4, c, ω_deg, "$(name): ω [°]",
                     output_dir, run_prefix, fig_post;
                     input_val=iv_ω, median_val=mv_ω)
        param_panel!(fig_post, row, 5, c, Ω_deg, "$(name): Ω [°]",
                     output_dir, run_prefix, fig_post;
                     input_val=iv_Ω, median_val=mv_Ω)    
    end
    
    output_path = joinpath(output_dir, "$(run_prefix)_posterior_batch.png")
    save(output_path, fig_post, px_per_unit=3)
    println("  Saved: $output_path")
    
    fig_post = nothing; GC.gc()
end

"""
Main function.

Processes all posterior runs in a directory and generates batch plots
comparing input parameters to recovered posterior distributions.
"""
function main(directory_path::String)
    isdir(directory_path) || error("Directory not found: $directory_path")
    
    
    all_files = readdir(directory_path, join=true)
    summaries = filter(f -> endswith(f, "_summary.md"), all_files)

    isempty(summaries) && error("No _summary.md files found in $directory_path")
    println("Found $(length(summaries)) run(s) to process.")
    
    processed = 0
    for summary_path in sort(summaries)
        run_dir = dirname(summary_path)
        run_prefix = replace(basename(summary_path), "_summary.md" => "")
        
        println("\nProcessing: $run_prefix")
        try
            summary_text = read(summary_path, String)
            cfg = extract_toml_from_summary(summary_text)
            
            # Exact epoch handling from plot_chain.jl
            epoch_mjd = get_epoch_mjd(cfg)
            epoch_year = cfg["epoch"]["year"]
            
            
            input_vals = get_input_parameters(cfg, epoch_mjd, stars)
            generate_posterior_plot(run_dir, run_prefix, input_vals, stars, output_dir)
            
            processed += 1
        catch e
            println("  *** ERROR: $e")
            showerror(stdout, e, backtrace())
        end
    end
    
    println("\nDone. Processed $processed/$(length(summaries)) runs.")
end

# ===================== RUN BATCH PLOTS =====================
if run_directory != ""
    main(run_directory)
elseif length(ARGS) >= 1
    main(ARGS[1])
else
    error("Usage: julia batch_posterior_plots.jl <directory>\nOR set run_directory variable.")
end