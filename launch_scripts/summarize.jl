"""
Summarize mock-data Octofitter results

Computes the median and 68% credible intervals overal from the median system parameters and semi-major axes for each run,
then compares the distribution of recovered posteriors with the inputted values.

This is used to assess the consistency of Ocen_IMBH_analysis/Octofitter pipeline 
and look for biases in parameter recovery. 

*** Note, this script parses the posterior_stats.txt files from each run to get the recovered posteriors, 
    it is not extracting results from the chain. 
"""

using TOML
using DataFrames
using Statistics
using CairoMakie

# ===================== CONFIGURATION =====================

# Select a directory with a bunch of runs to be summarized
# Make sure all the runs in this directory have the same config file, or face the error message. 
run_dir = raw"C:\Users\macke\Clusters\Ocen_IMBH_Mock_Analysis\mock_results\noOneil_16r_192c_June24";

outname = basename(run_dir)
outdir = joinpath(run_dir, outname)
mkpath(outdir)

# Stars
stars = ["A", "C", "D", "E", "F"]

# Colors
w = Makie.wong_colors()

# Find summary and posterior stats files for each run ID in run_dir
summary_files = filter(f -> endswith(f, "_summary.md"),
    readdir(run_dir, join=true))

posterior_files = filter(f -> endswith(f, "_posterior_stats.txt"),
    readdir(run_dir, join=true))

function runid(path)
    base = basename(path)
    replace(base,
        "_summary.md" => "",
        "_posterior_stats.txt" => "")
end

posterior_map = Dict(runid(f) => f for f in posterior_files)

# ===================== FUNCTIONS =====================

"""
Extracts the content of the first TOML code block from the provided string `text`. 
Used for getting the imput parameter configuration.
"""
function extract_toml_block(text::String)
    m = match(r"```toml(.*?)```"s, text)
    m === nothing && error("No TOML block found in summary file")
    return strip(m.captures[1])
end

"""
Parse a run summary file and extract the input mock parameters.
Returns a dictionary containing:
- System parameters: M_IMBH, plx, offsetx, offsety
- Derived semi-major axes for each star
"""
function parse_input_vals(summary_file)
    txt = read(summary_file, String)
    toml_txt = extract_toml_block(txt)
    cfg = TOML.parse(toml_txt)

    mock = cfg["mock"]
    
    # Extract system parameters
    M_IMBH = mock["M_IMBH"]
    
    # Calculate semi-major axes from orbital periods using Kepler's third law
    # a = cbrt(M * P^2) where P is period [yrs] and M is central mass [M_sun]
    A_period = mock["stars"]["A"]["orbital_elements"][1]
    C_period = mock["stars"]["C"]["orbital_elements"][1]  
    D_period = mock["stars"]["D"]["orbital_elements"][1]
    E_period = mock["stars"]["E"]["orbital_elements"][1]
    F_period = mock["stars"]["F"]["orbital_elements"][1]
    
    A_a = cbrt(M_IMBH * A_period^2)
    C_a = cbrt(M_IMBH * C_period^2)
    D_a = cbrt(M_IMBH * D_period^2)
    E_a = cbrt(M_IMBH * E_period^2)
    F_a = cbrt(M_IMBH * F_period^2)

    return Dict(
        "M_IMBH"  => M_IMBH,
        "plx"     => mock["plx"],
        "offsetx" => mock["offsetx"],
        "offsety" => mock["offsety"],
        # Calculated semi-major axes:
        "A_a" => A_a,
        "C_a" => C_a,
        "D_a" => D_a,
        "E_a" => E_a,
        "F_a" => F_a
    )
end

"""
Parse posterior summary statistics from a text file.
Extracts median and 68% CI of the posterior distribution for:
- M_IMBH [M_sun]
- plx, offsetx [mas], offsety [mas]
- Semi-major axes for each star [AU]
"""
function parse_posterior(file)

    lines = readlines(file)

    out = Dict{String, NamedTuple}()

    for line in lines
        line = strip(line)


        # IMBH mass 
        if startswith(line, "M_IMBH")

            m = match(r"M_IMBH.*?([\d\.\-eE]+)\s+\[\s*([\d\.\-eE]+),\s*([\d\.\-eE]+)\]",
                      line)

            if m !== nothing
                med = parse(Float64, m.captures[1]) * 1e4
                lo  = parse(Float64, m.captures[2]) * 1e4
                hi  = parse(Float64, m.captures[3]) * 1e4

                out["M_IMBH"] = (median=med, lo=lo, hi=hi)
            end
        end

        # Other system parameters
        for p in ["plx", "offsetx", "offsety"]

            startswith(line, p) 

            m = match(Regex("^$p.*?([\\d\\.\\-eE]+)\\s+\\[\\s*([\\d\\.\\-eE]+),\\s*([\\d\\.\\-eE]+)\\]"),
                      line)

            if m !== nothing
                out[p] = (
                    median=parse(Float64, m.captures[1]),
                    lo=parse(Float64, m.captures[2]),
                    hi=parse(Float64, m.captures[3])
                )
            end
        end        

        # Semi-major axes for each star 
        for star in stars
            pattern = Regex("^$(star): a \\[AU\\]\\s+([\\d\\.\\-eE]+)\\s+\\[\\s*([\\d\\.\\-eE]+),\\s*([\\d\\.\\-eE]+)\\]")
            m = match(pattern, line)
            if m !== nothing
                out["$(star)_a"] = (
                    median=parse(Float64, m.captures[1]),
                    lo=parse(Float64, m.captures[2]),
                    hi=parse(Float64, m.captures[3])
                )
            end
        end
    
    end
    return out
end

# ===================== LOAD INPUT PARAMETERS =====================

input_vals_list = Dict{String, Float64}[]


for sf in summary_files
    rid = runid(sf)
    push!(input_vals_list, parse_input_vals(sf))
end

if length(input_vals_list) == 0
    error("No valid runs found.")
end

input0 = input_vals_list[1]

for t in input_vals_list
    for k in keys(input0)
        if !isapprox(t[k], input0[k]; rtol=1e-12)
            error("Inconsistent input values across runs! Be more careful bro.")
        end
    end
end

# Determine stars from the first summary file
first_summary_file = summary_files[1]

println("Loaded $(length(summary_files)) runs")

# ===================== BUILD SUMMARY TABLE =====================

run_ids = String[]

M_true  = Float64[]
M_med   = Union{Float64,Missing}[]
M_lo    = Union{Float64,Missing}[]
M_hi    = Union{Float64,Missing}[]

plx_true = Float64[]
plx_med  = Union{Float64,Missing}[]
plx_lo   = Union{Float64,Missing}[]
plx_hi   = Union{Float64,Missing}[]

x_true = Float64[]
x_med  = Union{Float64,Missing}[]
x_lo   = Union{Float64,Missing}[]
x_hi   = Union{Float64,Missing}[]

y_true = Float64[]
y_med  = Union{Float64,Missing}[]
y_lo   = Union{Float64,Missing}[]
y_hi   = Union{Float64,Missing}[]

# Semi-major axis for each star
a_A_true = Float64[]
a_A_med  = Union{Float64,Missing}[]
a_A_lo   = Union{Float64,Missing}[]
a_A_hi   = Union{Float64,Missing}[]

a_C_true = Float64[]
a_C_med  = Union{Float64,Missing}[]
a_C_lo   = Union{Float64,Missing}[]
a_C_hi   = Union{Float64,Missing}[]

a_D_true = Float64[]
a_D_med  = Union{Float64,Missing}[]
a_D_lo   = Union{Float64,Missing}[]
a_D_hi   = Union{Float64,Missing}[]

a_E_true = Float64[]
a_E_med  = Union{Float64,Missing}[]
a_E_lo   = Union{Float64,Missing}[]
a_E_hi   = Union{Float64,Missing}[]

a_F_true = Float64[]
a_F_med  = Union{Float64,Missing}[]
a_F_lo   = Union{Float64,Missing}[]
a_F_hi   = Union{Float64,Missing}[]

# loop through summary files and posterior stats files
for sf in summary_files

    rid = runid(sf)

    input_vals = parse_input_vals(sf)
    post  = parse_posterior(posterior_map[rid])

    push!(run_ids, rid)

    # ---------------- M_IMBH ----------------
    push!(M_true, input_vals["M_IMBH"])

    if haskey(post, "M_IMBH")
        push!(M_med, post["M_IMBH"].median)
        push!(M_lo,  post["M_IMBH"].lo)
        push!(M_hi,  post["M_IMBH"].hi)
    else
        push!(M_med, missing); push!(M_lo, missing); push!(M_hi, missing)
    end

    # ---------------- plx ----------------
    push!(plx_true, input_vals["plx"])

    if haskey(post, "plx")
        push!(plx_med, post["plx"].median)
        push!(plx_lo,  post["plx"].lo)
        push!(plx_hi,  post["plx"].hi)
    else
        push!(plx_med, missing); push!(plx_lo, missing); push!(plx_hi, missing)
    end

    # ---------------- offsetx ----------------
    push!(x_true, input_vals["offsetx"])

    if haskey(post, "offsetx")
        push!(x_med, post["offsetx"].median)
        push!(x_lo,  post["offsetx"].lo)
        push!(x_hi,  post["offsetx"].hi)
    else
        push!(x_med, missing); push!(x_lo, missing); push!(x_hi, missing)
    end

    # ---------------- offsety ----------------
    push!(y_true, input_vals["offsety"])

    if haskey(post, "offsety")
        push!(y_med, post["offsety"].median)
        push!(y_lo,  post["offsety"].lo)
        push!(y_hi,  post["offsety"].hi)
    else
        push!(y_med, missing); push!(y_lo, missing); push!(y_hi, missing)
    end
    
    # ---------------- semi-major axis ----------------
    push!(a_A_true, input_vals["A_a"])
    if haskey(post, "A_a")
        push!(a_A_med, post["A_a"].median)
        push!(a_A_lo,  post["A_a"].lo)
        push!(a_A_hi,  post["A_a"].hi)
    else
        push!(a_A_med, missing); push!(a_A_lo, missing); push!(a_A_hi, missing)
    end

    push!(a_C_true, input_vals["C_a"])
    if haskey(post, "C_a")
        push!(a_C_med, post["C_a"].median)
        push!(a_C_lo,  post["C_a"].lo)
        push!(a_C_hi,  post["C_a"].hi)
    else
        push!(a_C_med, missing); push!(a_C_lo, missing); push!(a_C_hi, missing)
    end

    push!(a_D_true, input_vals["D_a"])
    if haskey(post, "D_a")
        push!(a_D_med, post["D_a"].median)
        push!(a_D_lo,  post["D_a"].lo)
        push!(a_D_hi,  post["D_a"].hi)
    else
        push!(a_D_med, missing); push!(a_D_lo, missing); push!(a_D_hi, missing)
    end

    push!(a_E_true, input_vals["E_a"])
    if haskey(post, "E_a")
        push!(a_E_med, post["E_a"].median)
        push!(a_E_lo,  post["E_a"].lo)
        push!(a_E_hi,  post["E_a"].hi)    
    else
        push!(a_E_med, missing); push!(a_E_lo, missing); push!(a_E_hi, missing)
    end

    push!(a_F_true, input_vals["F_a"])
    if haskey(post, "F_a")
        push!(a_F_med, post["F_a"].median)
        push!(a_F_lo,  post["F_a"].lo)
        push!(a_F_hi,  post["F_a"].hi)
    else
        push!(a_F_med, missing); push!(a_F_lo, missing); push!(a_F_hi, missing)
    end
end

rows = DataFrame(
    run_id = run_ids,

    M_true = M_true,
    M_med  = M_med,
    M_lo   = M_lo,
    M_hi   = M_hi,

    plx_true = plx_true,
    plx_med  = plx_med,
    plx_lo   = plx_lo,
    plx_hi   = plx_hi,

    offsetx_true = x_true,
    offsetx_med  = x_med,
    offsetx_lo   = x_lo,
    offsetx_hi   = x_hi,

    offsety_true = y_true,
    offsety_med  = y_med,
    offsety_lo   = y_lo,
    offsety_hi   = y_hi,
    
    A_a_true = a_A_true,
    A_a_med  = a_A_med,
    A_a_lo   = a_A_lo,
    A_a_hi   = a_A_hi,
    
    C_a_true = a_C_true,
    C_a_med  = a_C_med,
    C_a_lo   = a_C_lo,
    C_a_hi   = a_C_hi,
    
    D_a_true = a_D_true,
    D_a_med  = a_D_med,
    D_a_lo   = a_D_lo,
    D_a_hi   = a_D_hi,
    
    E_a_true = a_E_true,
    E_a_med  = a_E_med,
    E_a_lo   = a_E_lo,
    E_a_hi   = a_E_hi,
    
    F_a_true = a_F_true,
    F_a_med  = a_F_med,
    F_a_lo   = a_F_lo,
    F_a_hi   = a_F_hi
)

# ===================== SUMMARY STATISTICS =====================

summary = DataFrame(
    parameter = String[],
    input = Float64[],
    N = Int[],
    median_of_medians = Float64[],
    mean = Float64[],
    std = Float64[]
)

params = [
    ("M_IMBH", "M"),
    ("plx", "plx"),
    ("offsetx", "offsetx"),
    ("offsety", "offsety"),
    ("A_a", "A_a"),
    ("C_a", "C_a"),
    ("D_a", "D_a"),
    ("E_a", "E_a"),
    ("F_a", "F_a")
]

for (p, prefix) in params

    med = rows[!, Symbol(prefix * "_med")]
    lo  = rows[!, Symbol(prefix * "_lo")]
    hi  = rows[!, Symbol(prefix * "_hi")]

    input_val = input0[p]

    N = length(med)

    push!(summary, (
        p,
        input_val,
        N,
        median(med),
        mean(med),
        std(med),
    ))
end

# ===================== HISTOGRAMS =====================

plot_params = [
    ("M_IMBH", "M"),
    ("plx", "plx"),
    ("offsetx", "offsetx"),
    ("offsety", "offsety"),
    ("A_a", "A_a"),
    ("C_a", "C_a"),
    ("D_a", "D_a"),
    ("E_a", "E_a"),
    ("F_a", "F_a")
]

for (p, prefix) in plot_params

    med = rows[!, Symbol(prefix * "_med")]
    lo  = rows[!, Symbol(prefix * "_lo")]
    hi  = rows[!, Symbol(prefix * "_hi")]

    input = input0[p]

    # ---------------- histogram ----------------
    fig = Figure(size = (700, 450), figure_padding = 30)
    
    input_color = w[4]
    med_color  = w[3]
    hist_color = w[5]
    
    ax = Axis(fig[1, 1],
        xlabel=p,
        ylabel="count",
        title="Recovered $(p)"
    )

    hist!(ax, med; bins=25, color=(hist_color, 0.6), strokewidth = 0)


    m = median(med)

    vlines!(ax, [input, m],
        color=[input_color, med_color],
        linewidth=3
    )

    elements = [
        LineElement(color=input_color, linewidth=3),
        LineElement(color=med_color, linewidth=3),
    ]

    labels = ["Input value", "Median recovered"]

    axislegend(ax, elements, labels, position=:rb)
    ylims!(ax, 0, nothing)
    save(joinpath(outdir, "hist_$(p).png"), fig)

end

# ===================== SUMMARY MARKDOWN FILE =====================

report = joinpath(outdir, "$(outname).md")
open(report, "w") do io

    println(io, "# $outname \n")

    println(io, "Number of runs: $(nrow(rows))\n")

    println(io, "## Results Summary\n")

    # Separate system parameters from semi-major axes
    system_params = ["M_IMBH", "plx", "offsetx", "offsety"]
    sma_params = ["A_a", "C_a", "D_a", "E_a", "F_a"]

    # System parameters 
    println(io, "### System Parameters\n")
    for p in system_params
        r = summary[summary.parameter .== p, :][1, :]  
        println(io, "#### $(r.parameter)")
        println(io, "- Input value: $(input0[p])")
        println(io, "- Median recovered: $(r.median_of_medians)")
        println(io, "- Mean recovered: $(r.mean)")
        println(io, "- Std of recovered: $(r.std)")
        println(io, "")
    end

    # Semi-major axes 
    println(io, "### Semi-Major Axes\n")
    for p in sma_params
        r = summary[summary.parameter .== p, :][1, :]  
        println(io, "#### $(r.parameter)")
        println(io, "- Input value: $(input0[p])")
        println(io, "- Median recovered: $(r.median_of_medians)")
        println(io, "- Mean recovered: $(r.mean)")
        println(io, "- Std of recovered: $(r.std)")
        println(io, "")
    end
end

println("Outputs written to:")
println(outdir)