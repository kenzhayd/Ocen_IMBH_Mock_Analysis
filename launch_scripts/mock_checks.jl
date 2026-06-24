"""
This script checks orbital parameters being extracted from the real-data fitted chain.

Genereates orbits from orbital and system parameters used a TOML config file,
and compares the star'sorbit model position with the observed position. 

Useful for debugging when units, parameter names, coordinate system transformations get confusing.
"""

using Octofitter
using Octofitter: @variables
using CairoMakie
using Unitful
using UnitfulAstro
using TOML

include(joinpath(@__DIR__, "parse_config.jl"))
include("octo_utils.jl")

# ===================== CONFIGURATION =====================

# Load configuration
cfg = TOML.parsefile("C:\\Users\\macke\\Clusters\\Ocen_IMBH_Mock_Analysis\\configs\\mock_rand_3.toml")

# Load system parameters
epoch_year = cfg["epoch"]["year"]
epoch_mjd  = Octofitter.years2mjd(epoch_year)
mock_cfg = cfg["mock"]
M_IMBH  = mock_cfg["M_IMBH"]
plx     = mock_cfg["plx"]
offsetx = mock_cfg["offsetx"]
offsety = mock_cfg["offsety"]

# Colors
w = Makie.wong_colors()

# ===================== BIG ANALYSIS LOOP =====================
# Loop over all stars in config
for (name, star_cfg) in mock_cfg["stars"]

    println("Plotting star $name")

    # Extract orbital elements
    elems = star_cfg["orbital_elements"]
    P, e, i, ω, Ω, θ = elems

    # Load observed position from octo_utils.stars
    real = octo_utils.stars[name]
    obs_ra  = (real.ra  - octo_utils.ra_cm_deg)  * 3600 * 1000 * cosd(octo_utils.dec_cm_deg)
    obs_dec = (real.dec - octo_utils.dec_cm_deg) * 3600 * 1000


    # Build orbit from TOML file parameter set
    a_true = cbrt(M_IMBH * P^2)
    tp_true = θ_at_epoch_to_tperi(θ, epoch_mjd;
        a=a_true, e=e, i=i, ω=ω, Ω=Ω, M=M_IMBH)

    orbit_model = Visual{KepOrbit}(;
        a=a_true,
        e=e,
        i=i,
        ω=ω,
        Ω=Ω,
        tp=tp_true,
        M=M_IMBH,
        plx=plx
    )

    ts = range(epoch_mjd - P*365.25/2, epoch_mjd + P*365.25/2; length=600)

    ra_curve  = Float64[]
    dec_curve = Float64[]

    for t in ts
        sol = orbitsolve(orbit_model, t)
        ra  = raoff(sol) + offsetx
        dec = decoff(sol) + offsety
        push!(ra_curve, ra)
        push!(dec_curve, dec)
    end
    
    # Check if mock orbit passes through real position
    # Evaluate orbit model at the reference epoch
    sol_epoch = orbitsolve(orbit_model, epoch_mjd)
    ra_model_at_epoch  = raoff(sol_epoch) + offsetx
    dec_model_at_epoch = decoff(sol_epoch) + offsety

    # Compute residuals for positions at reference epoch (model  positions - observed position)
    d_ra  = ra_model_at_epoch  - obs_ra
    d_dec = dec_model_at_epoch - obs_dec

    println("\n--- Star $name diagnostics ---")
    println("Model orbit @ epoch:  dRA* = $(round(ra_model_at_epoch,  digits=4)) mas, dDec = $(round(dec_model_at_epoch, digits=4)) mas")
    println("Observed pos:          RA* = $(round(obs_ra,  digits=4)) mas, Dec = $(round(obs_dec, digits=4)) mas")
    println("Residual:             dRA* = $(round(d_ra, digits=4)) mas, dDec = $(round(d_dec, digits=4)) mas")

    # is the orbit centred on the IMBH correctly?
    ra_imbh_relative  = raoff(sol_epoch)   # star pos relative to IMBH only
    dec_imbh_relative = decoff(sol_epoch)  # star pos relative to IMBH only
    println("Star offset from IMBH: RA* offset = $(round(ra_imbh_relative, digits=4)) mas, Dec offset = $(round(dec_imbh_relative, digits=4)) mas")
    println("IMBH offset from cluster center: offsetx = $(round(offsetx, digits=8)) mas, offsety = $(round(offsety, digits=8)) mas")


    # Compare the real star's absolute RA/Dec to what stardata_struct would produce for mock positions
    real_ra_deg  = real.ra
    real_dec_deg = real.dec

    # What does stardata_struct compute for this mock star's absolute position?
    mas2deg = 1 / (3600 * 1000)
    mock_ra_deg  = octo_utils.ra_cm_deg  + (offsetx + ra_imbh_relative)  * mas2deg 
    mock_dec_deg = octo_utils.dec_cm_deg + (offsety + dec_imbh_relative) * mas2deg

    println("Real star RA/Dec (deg):  $(real_ra_deg),  $(real_dec_deg)")
    println("Mock star RA/Dec (deg):  $(round(mock_ra_deg, digits=8)),  $(round(mock_dec_deg, digits=8))")

    # Flag errors
    println("\nDo we still have errors?\n")
    if abs(d_ra) < 1.0 && abs(d_dec) < 1.0
        println("Mock orbit passes through observed position within 1 mas so maybe life's good?\n")
    elseif abs(d_ra - offsetx) < 1.0 || abs(d_dec - offsety) < 1.0
        println("Offset night not be applied: orbit could be centred on cluster center, not IMBH.\n")
    elseif abs(mock_ra_deg - real_ra_deg)*3600*1000 > 5.0
        println(" Mock and real star positions differ by >5 mas: something is fishy.\n")    
    else
        println("Maybe life's good?\n")
    end

    # ===================== PLOT CHECKS =====================

    fig = Figure(size=(650, 650))
    ax = Axis(fig[1,1];
        xlabel="dα* [mas]",
        ylabel="dd [mas]",
        aspect=DataAspect(),
        xreversed=true,
        title="Star $name"
    )

    CairoMakie.lines!(ax, ra_curve, dec_curve; color=:blue, linewidth=2, label="Mock Orbit")
    CairoMakie.scatter!(ax, [obs_ra], [obs_dec];
        color=:red, markersize=14, marker=:star5, label="Real Observed Position")

    CairoMakie.scatter!(ax, [offsetx], [offsety];
        color=:black, markersize=12, marker=:circle, label="IMBH")

    axislegend(ax; position=:rb)
    display(fig)
end




