"""
Extract samples from real-data fitted chain.

Prints TOML blocks of orbital parameters for a given sample,
used for mock recovery analysis config files.
"""

using Statistics
using Octofitter
using Printf
using Random

# =====================  CONFIGURATION ===================== 

chain_path = raw"C:\Users\macke\Clusters\run_108836842\starsACDEF_192c_18r_cont_10836842_chain.fits"
chain = Octofitter.loadchain(chain_path)
println("Loaded chain: $chain_path")

summary_path = raw"C:\Users\macke\Clusters\run_108836842\starsACDEF_192c_18r_cont_10836842_summary.md"
summary_text = read(summary_path, String)

println("\n=== ALL COLUMN NAMES IN CHAIN ===")
for c in sort(names(chain))
    println(c)
end
println("=================================\n")

# Star names
line = first(l for l in split(summary_text, '\n') if occursin("**Stars:**", l))
after_stars = split(line, "**Stars:**")[2]
star_names = strip.(split(strip(after_stars), ","))
println("Extracted star names: ", star_names)

# ===================== EXTRACT PARAMETERS ===================== 

col(chain, sym) = vec(Array(chain[sym]))

function get_sample_value(chain, idx, star, param)
    return col(chain, Symbol("$(star)_$(param)"))[idx]
end

function get_system_value(chain, idx, param)
    return col(chain, Symbol(param))[idx]
end


# ===================== PRINT CONFIG BLOCK ===================== 

function print_mock_toml_block(chain, star_names, idx;
                               target_M=nothing,
                               label="sample")

    Mval       = get_system_value(chain, idx, :M)
    plxval     = get_system_value(chain, idx, :plx)
    offsetxval = get_system_value(chain, idx, :offsetx)
    offsetyval = get_system_value(chain, idx, :offsety)

    println()
    println("# ========== MOCK DATA CONFIGURATION ==========")
    println()
    println("[mock]")
    println("enabled = true")
    println()

    println("# Mock central IMBH parameters")
    @printf("M_IMBH  = %.12f    # Solar masses\n", Mval)
    @printf("plx     = %.12f    # parallax [mas]\n", plxval)
    @printf("offsetx = %.12f\n", offsetxval)
    @printf("offsety = %.12f\n", offsetyval)
    println()

    println("$label")
    println("# Fitting should recover these values")
    println("# Sample index: $idx")

    if target_M !== nothing
        @printf("# Target central mass: M ≈ %.6f\n", target_M)
    end

    println("# period, eccentricity, inclination [rad], ω [rad], Ω [rad], theta [rad]")
    println()
    println("[mock.stars]")
    println()

    for star in star_names
        P = get_sample_value(chain, idx, star, "P")
        e = get_sample_value(chain, idx, star, "e")
        i = get_sample_value(chain, idx, star, "i")
        ω = get_sample_value(chain, idx, star, "ω")
        Ω = get_sample_value(chain, idx, star, "Ω")
        θ = get_sample_value(chain, idx, star, "θ")

        println("[mock.stars.$star]")
        @printf("orbital_elements = [%.12f, %.12f, %.12f, %.12f, %.12f, %.12f]\n",
                P, e, i, ω, Ω, θ)
        println()
    end
end

# ===================== RANDOM SAMPLES ===================== 

function random_mass_indices(Mvals, target_M; tol=600.0, n=3, rng=Random.default_rng())

    mask = abs.(Mvals .- target_M) .< tol
    idx = findall(mask)

    isempty(idx) && error("No samples found in mass window")

    idx = Random.shuffle(rng, idx)
    k = min(n, length(idx))
    return idx[1:k]
end

# ===================== EXTRACT SAMPLES FROM CHAIN ===================== 

logpost = col(chain, :logpost)
Mvals   = col(chain, :M)

orbit_params = ["P","a","e","i","ω","Ω","θ","tp"]

# ===================== MAP ===================== 

idx_MAP = argmax(logpost)
println("MAP sample index: $idx_MAP  |  log-posterior: $(logpost[idx_MAP])")

println("\n=== MAP orbital elements ===")
for star in star_names
    println("Star $star")
    println("-------------------------")
    for p in orbit_params
        try
            val = get_sample_value(chain, idx_MAP, star, p)
            println(rpad(p, 5), " = ", val)
        catch
        end
    end
    println()
end

println("=== MAP system parameters ===")
for p in ["M","plx","offsetx","offsety"]
    sym = Symbol(p)
    println(rpad(p, 8), " = ", col(chain, sym)[idx_MAP])
end

print_mock_toml_block(chain, star_names, idx_MAP;
                      label="MAP posterior sample")

# ===================== CLOSEST-MASS ===================== 

target_M = 64000.0
idx_M = argmin(abs.(Mvals .- target_M))

println("\nClosest sample for M=$(target_M):")
println("  Index = $idx_M,  M = $(Mvals[idx_M]),  ΔM = $(abs(Mvals[idx_M] - target_M))")

println("\n=== Closest orbital elements ===")
for star in star_names
    println("Star $star")
    println("-------------------------")
    for p in orbit_params
        try
            val = get_sample_value(chain, idx_M, star, p)
            println(rpad(p, 5), " = ", val)
        catch
        end
    end
    println()
end

println("=== Closest-mass system parameters ===")
for p in ["M","plx","offsetx","offsety"]
    sym = Symbol(p)
    println(rpad(p, 8), " = ", col(chain, sym)[idx_M])
end

print_mock_toml_block(chain, star_names, idx_M;
                      target_M=target_M,
                      label="Closest mass sample")


# ===================== MAX LOGPOST, MASS WITHIN RANGE ===================== 

tol = 600.0
mask = abs.(Mvals .- target_M) .< tol
idx_candidates = findall(mask)

isempty(idx_candidates) && error("No candidates in mass window")

idx_best_logpost = idx_candidates[argmax(logpost[idx_candidates])]

println("\nMax logpost sample for M=$(target_M):")
println("  Index = $idx_best_logpost,  M = $(Mvals[idx_best_logpost]),  ΔM = $(abs(Mvals[idx_best_logpost] - target_M))")

println("\n=== Max logpost orbital elements ===")
for star in star_names
    println("Star $star")
    println("-------------------------")
    for p in orbit_params
        try
            val = get_sample_value(chain, idx_best_logpost, star, p)
            println(rpad(p, 5), " = ", val)
        catch
        end
    end
    println()
end

println("=== Max logpost system parameters ===")
for p in ["M","plx","offsetx","offsety"]
    sym = Symbol(p)
    println(rpad(p, 8), " = ", col(chain, sym)[idx_best_logpost])
end

print_mock_toml_block(chain, star_names, idx_best_logpost;
                      target_M=target_M,
                      label="Max logpost mass sample")

# ===================== UNIFORMLY SAMPLED INDICES ===================== 

tol_random = 600.0
n_random = 3

idx_random = random_mass_indices(Mvals, target_M;
                                 tol=tol_random,
                                 n=n_random)

println("\n=== Random samples near M=$(target_M) (±$(tol_random) Msun) ===")

for (k, idx) in enumerate(idx_random)

    M_sample = Mvals[idx]
    ΔM = abs(M_sample - target_M)

    println("\n--- Random sample $(k) ---")
    println("Index = $idx,  M = $M_sample,  ΔM = $ΔM")

    println("\nOrbital elements for random sample $(k):")
    for star in star_names
        println("Star $star")
        println("-------------------------")
        for p in orbit_params
            try
                val = get_sample_value(chain, idx, star, p)
                println(rpad(p, 5), " = ", val)
            catch
            end
        end
        println()
    end

    println("System parameters for random sample $(k):")
    for p in ["M","plx","offsetx","offsety"]
        sym = Symbol(p)
        println(rpad(p, 8), " = ", col(chain, sym)[idx])
    end

    label = @sprintf("Random mass sample %d (M ≈ %.3f, ΔM ≈ %.3f)",
                     k, M_sample, ΔM)

    print_mock_toml_block(chain, star_names, idx;
                          target_M=target_M,
                          label=label)
end