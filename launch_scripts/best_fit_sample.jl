using Statistics
using Octofitter
using Printf

# ── 3. Load chain ────────────────────────────────────────────────────────────
chain_path = "C:\\Users\\macke\\Clusters\\Ocen_IMBH_Mock_Analysis\\108836842\\starsACDEF_192c_18r_cont_10836842_chain.fits"

chain = Octofitter.loadchain(chain_path)
println("Loaded chain: $chain_path")
println(chain)
stats_path = joinpath("C:\\Users\\macke\\Clusters\\Ocen_IMBH_Mock_Analysis\\108836842\\starsACDEF_192c_18r_cont_10836842_posterior_stats.txt")

# ── 4. Get star names from summary, verify against chain columns ──────────────
# FITS column names are ASCII-only; Unicode characters (ω, Ω) may be encoded
# differently on save/load.  Use the summary as the authoritative source for
# star names, then discover the actual column names for each orbital element.

col_names = Set(Symbol.(names(chain)))
summary_text = read("C:\\Users\\macke\\Clusters\\Ocen_IMBH_Mock_Analysis\\108836842\\starsACDEF_192c_18r_cont_10836842_summary.md", String)
star_names = sort!(String.(strip.(split(summary_stars_line[1], ","))))



orbit_params = ["a","e","i","ω","Ω","tp"]




# Good mass sample

target_M = 64000.0
Mvals = vec(Array(chain[:M]))

idx = argmin(abs.(Mvals .- target_M))

println("Closest sample:")
println("Index = ", idx)
println("M = ", Mvals[idx])
println("ΔM = ", abs(Mvals[idx] - target_M))

for star in star_names
    println("Star $star")
    println("-------------------------")

    for p in orbit_params
        try
            col = find_col(col_names, star, p)
            val = chain[col][idx]
            println(rpad(p,5), " = ", val)
        catch
            # skip missing params
        end
    end

    println()
end


# === System parameters ===

param_names = Set(Symbol.(names(chain)))

println("System Parameters")
println("-------------------------")

for p in system_params
    sym = Symbol(p)
    if sym in param_names
        vals = vec(Array(chain[sym]))
        println(rpad(p,8), " = ", vals[idx])
    else
        println(rpad(p,8), " = MISSING")
    end
end
8




# Map sample
logpost = vec(Array(chain[:logpost]))
idx_MAP = argmax(logpost)
for star in star_names
    println("Star $star")
    println("-------------------------")

    for p in orbit_params
        try
            col = find_col(col_names, star, p)
            val = chain[col][idx_MAP]
            println(rpad(p,5), " = ", val)
        catch
            # skip missing params
        end
    end

    println()
end


# === System parameters ===

param_names = Set(Symbol.(names(chain)))

println("System Parameters")
println("-------------------------")

for p in system_params
    sym = Symbol(p)
    if sym in param_names
        vals = vec(Array(chain[sym]))
        println(rpad(p,8), " = ", vals[idx_MAP])
    else
        println(rpad(p,8), " = MISSING")
    end
end