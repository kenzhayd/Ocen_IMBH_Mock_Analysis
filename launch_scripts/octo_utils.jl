module octo_utils

# ========== Environment variables ==========
ENV["JULIA_NUM_THREADS"] = "auto"
ENV["OCTOFITTERPY_AUTOLOAD_EXTENSIONS"] = "yes"

# ========== Imports ==========

using Unitful
using UnitfulAstro
using LinearAlgebra
using Octofitter
using OctofitterRadialVelocity
using Distributions
using Statistics



# ========== Define StarData Type ==========
"""
Utilities for organizing and fitting multi-star astrometric data with Octofitter.
Includes:
- StarData struct for storing motion and error data
- stars dictionary with input for each star
- Functions to build models, propagate errors, and run Octofitter fits and plots
"""
struct StarData
    name::String
    ra::Float64         # Right Ascension position in degrees
    dec::Float64        # Declination position in degrees
    sigma_ra::Float64   # RA position uncertainty in mas
    sigma_dec::Float64  # Dec position uncertainty in mas
    pm_ra::Float64      # Proper motion in RA (mas/year)
    pm_dec::Float64     # Proper motion in Dec (mas/year)
    acc_ra::Float64     # Acceleration in RA (mas/year²)
    acc_dec::Float64    # Acceleration in Dec (mas/year²)
    sigma_pm_ra::Float64
    sigma_pm_dec::Float64
    sigma_acc_ra::Float64
    sigma_acc_dec::Float64
    v2D::Float64
    v2D_err::Float64
    rv::Float64
    rv_err::Float64
end

# ========== Omega Centauri Data ==========
# Cluster center from Anderson & van der Marel (2010), ApJ 710, 1032:
# (α, δ) = (13:26:47.24, −47:28:46.45)  →  201.6968333°, −47.4795694°
# https://iopscience.iop.org/article/10.1088/0004-637X/710/2/1032
# Used as the origin of the relative astrometry frame; the IMBH position
# is a free parameter (offsetx, offsety) relative to this point.

# Center of mass RA and Dec in deg
ra_cm_deg  = 201.6968333
dec_cm_deg = -47.4795694

# Distance to Omega Centauri center (kiloparsecs and converted to km)
distance_kpc = 5.43u"kpc"
distance_km = uconvert(u"km", distance_kpc)

# Assumed errors in position (mas)
ra_err = 0.5u"mas"
dec_err = 0.5u"mas"

# Cluster systemic radial velocity (Baumgardt catalogue)
const rv_cluster     = 232780.0   # [m/s]  (232.78 ± 0.21 km/s)
const rv_cluster_err = 210.0      # [m/s]

# Cluster central escape velocity (Häberle et al. 2024)
const v_esc_cluster_kms = 62.0    # [km/s]

# ========== Define Error Propagation Function ==========
"""
Propagates uncertainty in position due to uncertainties in
proper motion and acceleration over time. No uncertainty in initial position is available.

Parameters:
- sigma_pos : Initial uncertainty in position (mas)
- sigma_pm  : Uncertainty in proper motion (mas/yr)
- sigma_acc : Uncertainty in acceleration (mas/yr²)
- dt        : Time from reference epoch (in years)

Returns:
- sigma_pos : Total uncertainty in predicted position at time dt (mas)
"""
function propagate_error(sigma_pos, sigma_pm, sigma_acc, dt)
    term_pos = sigma_pos^2                      
    term_pm = (dt * sigma_pm)^2                 
    term_acc = (0.5 * dt^2 * sigma_acc)^2      

    sqrt(term_pos + term_pm + term_acc)
end

# ========== Define Position Projection Function ==========
"""
Gives a fake observed position using an observed angular position, velocity, and acceleration

Parameters:
- pos : Initial position (in mas)
- pm : Proper motion (in mas/yr)
- acc : Acceleration (in mas/yr²)
- dt : Time offset(s) from the reference epoch in years

Returns:
- pos_final : calculated position (in mas)
"""
function fake_pos(pos, pm, acc, dt)
    pos + pm*dt + 0.5*acc*dt^2
end

# ========== Define Star Dictionary ==========
stars = Dict{String,StarData}(
    "A" => StarData(
        "A",
        201.6967263, # deg
        -47.4795835, # deg  
        ustrip(ra_err), # mas
        ustrip(dec_err), # mas
        3.563, # mas/year
        2.564, # mas/year
        -0.0069, # mas/year²
        0.0085,  # mas/year²
        0.038,   # mas/year
        0.055,   # mas/year
        0.0083,  # mas/year²
        0.0098,  # mas/year²
        113.0,   # km/s
        1.1,      # km/s error
        NaN,
        NaN
    ),
    "B" => StarData(
        "B",
        201.6968888,
        -47.4797138,
        ustrip(ra_err), # mas
        ustrip(dec_err), # mas
        2.167,
        1.415,
        0.0702,
        0.0228,
        0.182,
        0.081,
        0.0239,
        0.0157,
        66.6,
        4.1,
        NaN,
        NaN
    ),
    "C" => StarData(
        "C",
        201.6966378,
        -47.4793672,
        ustrip(ra_err), # mas
        ustrip(dec_err), # mas
        1.117,
        3.514,
        0.0028,
        -0.0060,
        0.127,
        0.056,
        0.0333,
        0.0123,
        94.9,
        1.7,
        NaN,
        NaN
    ),
    "D" => StarData(
        "D",
        201.6968346,
        -47.4793233,
        ustrip(ra_err), # mas
        ustrip(dec_err), # mas
        2.559,
        -1.617,
        0.0357,
        -0.0194,
        0.082,
        0.061,
        0.0177,
        0.0162,
        77.9,
        2.0,
        NaN,
        NaN
    ),
    "E" => StarData(
        "E",
        201.6973080,
        -47.4797545,
        ustrip(ra_err), # mas
        ustrip(dec_err), # mas
        -2.149,
        1.638,
        0.0072,
        -0.0009,
        0.025,
        0.037,
        0.0042,
        0.0075,
        69.6,
        0.8,
        261700.0,
        2700.0
    ),
    "F" => StarData(
        "F",
        201.6977125,
        -47.4792625,
        ustrip(ra_err), # mas
        ustrip(dec_err), # mas
        0.436,
        -2.584,
        0.0052,
        -0.0015,
        0.017,
        0.016,
        0.0038,
        0.0038,
        67.4,
        0.4,
        232500.0,
        4000.0
        ),
    "G" => StarData(
        "G",
        201.6961340,
        -47.4790585,
        ustrip(ra_err), # mas
        ustrip(dec_err), # mas
        -1.317,
        2.207,
        -0.0197,
        0.0173,
        0.098,
        0.062,
        0.0267,
        0.0170,
        66.2,
        1.9,
        NaN,
        NaN
    ),
)

# ========== Astrometry Input for Octofitter ==========
"""
Simulates astrometric positions for a star at three epochs: past, present, and future.

Parameters:
- star: StarData object containing motion and error data
- epoch: Central epoch in calendar years (e.g., 2010)
- dt: Time offset from the central epoch in years (e.g., 10)

Returns:
- epochs_mjd: Modified Julian Dates for the three epochs
- ra_rel: Relative right ascension values at the three epochs (mas)
- dec_rel: Relative declination values at the three epochs (mas)
- ra_errs: Measurement errors for RA at the three epochs (mas)
- dec_errs: Measurement errors for Dec at the three epochs (mas)
"""

function simulate_astrometry(star::StarData, epoch::Real, dt::Real)
    # Epochs in years and MJD
    epochs_years = [epoch - dt, epoch, epoch + dt]
    epochs_mjd = [Octofitter.years2mjd(y) for y in epochs_years]

    # --- Absolute positions in mas ---
    ra0_mas  = star.ra * 3600 * 1000    
    dec0_mas = star.dec * 3600 * 1000

    # --- Propagate positions using proper motion & acceleration  ---
    past_ra_mas   = fake_pos(ra0_mas, star.pm_ra, star.acc_ra, -dt)
    past_dec_mas  = fake_pos(dec0_mas, star.pm_dec, star.acc_dec, -dt)
    future_ra_mas = fake_pos(ra0_mas, star.pm_ra, star.acc_ra, dt)
    future_dec_mas= fake_pos(dec0_mas, star.pm_dec, star.acc_dec, dt)

    ra_abs = [past_ra_mas, ra0_mas, future_ra_mas]
    dec_abs = [past_dec_mas, dec0_mas, future_dec_mas]

    # --- Convert to relative RA/Dec (Δα*, Δδ) ---
    # Apply cos(δ_ref) to RA differences to match the α* convention used by Octofitter.
    ra_rel  = (ra_abs .- (ra_cm_deg  * 3600 * 1000)) .* cosd(dec_cm_deg)
    dec_rel =  dec_abs .- (dec_cm_deg * 3600 * 1000)

    # --- Error propagation ---
    past_ra_err   = propagate_error(ustrip(ra_err), star.sigma_pm_ra, star.sigma_acc_ra, -dt)
    past_dec_err  = propagate_error(ustrip(dec_err), star.sigma_pm_dec, star.sigma_acc_dec, -dt)
    future_ra_err = propagate_error(ustrip(ra_err), star.sigma_pm_ra, star.sigma_acc_ra, dt)
    future_dec_err= propagate_error(ustrip(dec_err), star.sigma_pm_dec, star.sigma_acc_dec, dt)

    ra_errs  = [past_ra_err, ustrip(ra_err), future_ra_err]
    dec_errs = [past_dec_err, ustrip(dec_err), future_dec_err]

    return epochs_mjd, ra_rel, dec_rel, ra_errs, dec_errs
end






# ========================================================
#  Angular to Linear Acceleration Conversion
# ========================================================
"""
Convert angular acceleration from milliarcseconds per year squared (mas/yr²)
to linear acceleration in kilometers per second squared (km/s²).

# Arguments
- a_masyr2::Quantity : Angular acceleration with units mas/yr²
- distance_km::Quantity : Distance to the object with units km

# Returns
- Linear acceleration in km/s² as a Quantity
"""
function masyr2_to_kms2(a_masyr2::Unitful.Quantity, distance_km::Unitful.Quantity)
    # 1 mas = 1e-3 arcsec = (1e-3 / 3600) deg = (1e-3 / 3600)*(π/180) rad
    # convert manually:
    a_radyr2 = uconvert(u"rad"/u"yr"^2, a_masyr2 * (1e-3 / 3600) * (π / 180))

    # Linear acceleration = angular acceleration × distance
    # a_radyr2 [rad/yr²] * distance_km [km] = km/yr²
    a_kmyr2 = a_radyr2 * distance_km

    # Convert km/yr² to km/s²
    a_kms2 = uconvert(u"km"/u"s"^2, a_kmyr2)

    return a_kms2
end

# ========================================================
#  Total Angular Accelerations and Uncertainty 
# ========================================================
"""
Compute total plane-of-sky angular and physical acceleration
(and their uncertainties) for a given star.

Parameters
----------
star : StarData
    A star object with acc_ra, acc_dec, sigma_acc_ra, and sigma_acc_dec attributes.

Returns
-------
a_total_masyr2 : Float64
    Total angular acceleration in mas/yr².
a_total_masyr2_err : Float64
    Uncertainty in total angular acceleration.
a_total_kms2 : Float64
    Total physical acceleration in km/s².
a_total_kms2_err : Float64
    Uncertainty in physical acceleration.
"""
function total_accelerations(star::StarData)
    a_ra = star.acc_ra
    a_dec = star.acc_dec
    a_ra_err = star.sigma_acc_ra
    a_dec_err = star.sigma_acc_dec

    # Angular acceleration magnitude
    a_total_masyr2 = sqrt(a_ra^2 + a_dec^2)

    # Uncertainty propagation
    a_total_masyr2_err = sqrt(
        (a_ra * a_ra_err / a_total_masyr2)^2 +
        (a_dec * a_dec_err / a_total_masyr2)^2
    )

    # Convert to physical acceleration [km/s²]
    a_total_kms2 = masyr2_to_kms2(a_total_masyr2 * u"mas/yr^2", distance_km)
    a_total_kms2_err = masyr2_to_kms2(a_total_masyr2_err * u"mas/yr^2", distance_km)

    return a_total_masyr2, a_total_masyr2_err, a_total_kms2, a_total_kms2_err
end

# ========================================================
#  Build Direct Observation Objects for Octofitter v8
# ========================================================
"""
Build observation objects from a StarData object at a single epoch.
Uses the direct kinematic observables rather than synthetic multi-epoch astrometry.

Parameters:
- star: StarData object with position, PM, acceleration, and (optional) RV data
- epoch_mjd: Observation epoch in Modified Julian Date
- include_rv: Whether to build RV observation (default true; only used if star has RV data)
- z_prior_sigma: If not nothing, build a PlanetZPriorObs with Normal(0, z_prior_sigma) in AU
- include_esc_vel: If true, build a PlanetEscVelObs using star.v2D / star.v2D_err [km/s]

Returns:
- astrom: PlanetRelAstromObs (single-epoch relative position)
- pm: PlanetPMObs (single-epoch proper motion)
- acc: PlanetAccelObs (single-epoch acceleration)
- rv: PlanetRelativeRVObs or nothing (single-epoch peculiar radial velocity)
- zp: PlanetZPriorObs or nothing (LOS position prior)
- ev: PlanetEscVelObs or nothing (Häberle-style escape velocity constraint)
"""
function build_star_observations(star::StarData, epoch_mjd::Float64;
                                  include_rv::Bool=true,
                                  z_prior_sigma::Union{Nothing,Float64}=nothing,
                                  include_esc_vel::Bool=false,
                                  acceleration_type::String="vector")
    # 1. Single-epoch position relative to cluster center.
    # RA offset is multiplied by cos(δ_ref) to give Δα* (east in mas), consistent
    # with the α* convention used by raoff(sol), pmra(sol), and the input PM/accel data.
    ra_rel_mas  = (star.ra  - ra_cm_deg)  * 3600 * 1000 * cosd(dec_cm_deg)
    dec_rel_mas = (star.dec - dec_cm_deg) * 3600 * 1000
    astrom = PlanetRelAstromObs(
        (epoch=[epoch_mjd], ra=[ra_rel_mas], dec=[dec_rel_mas],
         σ_ra=[ustrip(ra_err)], σ_dec=[ustrip(dec_err)], cor=[0.0]);
        name="$(star.name)_pos"
    )

    # 2. Proper motion at same epoch
    pm = PlanetPMObs(
        (epoch=[epoch_mjd], pmra=[star.pm_ra], pmdec=[star.pm_dec],
         σ_pmra=[star.sigma_pm_ra], σ_pmdec=[star.sigma_pm_dec], cor=[0.0]);
        name="$(star.name)_pm"
    )

    # 3. Acceleration at same epoch — vector components or magnitude only
    acc = if acceleration_type == "vector"
        PlanetAccelObs(
            (epoch=[epoch_mjd], accra=[star.acc_ra], accdec=[star.acc_dec],
             σ_accra=[star.sigma_acc_ra], σ_accdec=[star.sigma_acc_dec], cor=[0.0]);
            name="$(star.name)_acc"
        )
    elseif acceleration_type == "magnitude"
        # Propagate uncertainty: σ_|a| = √((aα·σ_aα)² + (aδ·σ_aδ)²) / |a|
        # Guard against zero magnitude (undetected acceleration): fall back to RSS of σ components.
        _accmag = hypot(star.acc_ra, star.acc_dec)
        _σ_accmag = _accmag > 0 ?
            sqrt((star.acc_ra * star.sigma_acc_ra)^2 + (star.acc_dec * star.sigma_acc_dec)^2) / _accmag :
            hypot(star.sigma_acc_ra, star.sigma_acc_dec)
        PlanetAccelMagObs(
            (epoch=[epoch_mjd], accmag=[_accmag], σ_accmag=[_σ_accmag]);
            name="$(star.name)_accmag"
        )
    else
        error("Unknown acceleration_type: \"$(acceleration_type)\". Expected \"vector\" or \"magnitude\".")
    end

    # 4. Radial velocity (peculiar, relative to cluster systemic RV)
    rv = nothing
    if include_rv && !isnan(star.rv) && !isnan(star.rv_err)
        rv_peculiar = star.rv - rv_cluster
        σ_rv_total  = hypot(star.rv_err, rv_cluster_err)
        rv = PlanetRelativeRVObs(
            (epoch = epoch_mjd, rv = rv_peculiar, σ_rv = σ_rv_total);
            name = "$(star.name)_rv",
            variables = @variables begin end
        )
    end

    # 5. LOS position prior (z ~ Normal(0, σ_z) in AU)
    zp = nothing
    if z_prior_sigma !== nothing
        zp = PlanetZPriorObs(epoch_mjd, Normal(0.0, z_prior_sigma);
                              name="$(star.name)_zprior")
    end

    # 6. Escape velocity constraint (Häberle et al. 2024)
    ev = nothing
    if include_esc_vel
        ev = PlanetEscVelObs(epoch_mjd, star.v2D, star.v2D_err, v_esc_cluster_kms;
                              name="$(star.name)_escvel")
    end

    return astrom, pm, acc, rv, zp, ev
end



# ========== START MOCK MODIFICATIONS  ==========
# The following functions have been added for generating mock data to test Octofitter fitting.
# See MOCK_PLAN.md for details.

# ========================================================
#  Mock Orbit Generation
# ========================================================
"""
Container for a mock orbit generated around a central mass.
"""
struct MockOrbit
    name::String
    orbit::Any  # Visual{KepOrbit} object
end

"""
Generate a Keplerian orbit for a mock star around a central mass.

Parameters:
- name: Star identifier 
- a: Semi-major axis [AU]
- e: Eccentricity
- i: Inclination [rad]
- ω: Argument of periastron [rad]
- Ω: Longitude of ascending node [rad]
- M: Central mass [solar masses]
- plx: Parallax [mas]
- t_ref: Reference epoch (MJD) for setting periastron time

Returns:
- MockOrbit containing a Visual{KepOrbit} object
"""
function make_star(name; a, e, i, ω, Ω, M, plx, t_ref)
    θ = 2π * rand()
    tp = θ_at_epoch_to_tperi(
        θ, t_ref;
        a=a, e=e, i=i, ω=ω, Ω=Ω, M=M
    )
    orbit = Visual{KepOrbit}(;
        a=a,
        e=e,
        i=i,
        ω=ω,
        Ω=Ω,
        M=M,
        tp=tp,
        plx=plx
    )
    return MockOrbit(name, orbit)
end

"""
Calculate noiseless data (position, proper motion, acceleration, RV) for a mock star.

Parameters:
- star: MockOrbit containing the Keplerian orbit
- epoch: Observation epoch (MJD)

Returns:
- NamedTuple with fields:
  - raoff: RA offset [mas]
  - decoff: Dec offset [mas]
  - pmra: Proper motion in RA [mas/yr]
  - pmdec: Proper motion in Dec [mas/yr]
  - accra: Acceleration in RA [mas/yr²]
  - accdec: Acceleration in Dec [mas/yr²]
  - rv: Line-of-sight radial velocity [m/s]
"""
function mock_data(star::MockOrbit, epoch::Real)
    sol = orbitsolve(star.orbit, epoch)
    return (
        raoff     = raoff(sol),
        decoff    = decoff(sol),
        pmra      = pmra(sol),
        pmdec     = pmdec(sol),
        accra     = accra(sol),
        accdec    = accdec(sol),
        rv        = radvel(sol)
    )
end

"""
Generate a mock StarData object from orbital parameters.

Creates a Keplerian orbit, evaluates noiseless observables at a reference epoch,
and packages them into a StarData struct with specified observation uncertainties.

Parameters:
- name: Star name/identifier
- a: Semi-major axis [AU]
- e: Eccentricity
- i: Inclination [rad]
- ω: Argument of periastron [rad]
- Ω: Longitude of ascending node [rad]
- M: Central mass [solar masses]
- plx: Parallax [mas]
- t_ref: Reference epoch for orbit generation (MJD)
- epoch: Epoch at which to evaluate observables (MJD)
- sigma_ra: Uncertainty in RA position [mas]
- sigma_dec: Uncertainty in Dec position [mas]
- sigma_pm_ra: Uncertainty in proper motion RA [mas/yr]
- sigma_pm_dec: Uncertainty in proper motion Dec [mas/yr]
- sigma_acc_ra: Uncertainty in acceleration RA [mas/yr²]
- sigma_acc_dec: Uncertainty in acceleration Dec [mas/yr²]
- sigma_rv: Uncertainty in radial velocity [m/s]

Returns:
- StarData struct with noiseless observables and specified uncertainties
"""
function stardata_struct(name;
    a, e, i, ω, Ω,
    M,
    plx,
    t_ref,
    epoch,
    sigma_ra,
    sigma_dec,
    sigma_pm_ra,
    sigma_pm_dec,
    sigma_acc_ra,
    sigma_acc_dec,
    sigma_rv
)
    orbit = make_star(name; a, e, i, ω, Ω, M=M, plx=plx, t_ref=t_ref)
    obs = mock_data(orbit, epoch)
    
    v2D = hypot(obs.pmra, obs.pmdec)

    v2D_err = sqrt(
        (obs.pmra / v2D)^2 * sigma_pm_ra^2 +
        (obs.pmdec / v2D)^2 * sigma_pm_dec^2
        )

    return StarData(
    name,
    obs.raoff,
    obs.decoff,
    sigma_ra,
    sigma_dec,
    obs.pmra,
    obs.pmdec,
    obs.accra,
    obs.accdec,
    sigma_pm_ra,
    sigma_pm_dec,
    sigma_acc_ra,
    sigma_acc_dec,
    v2D,
    v2D_err,
    obs.rv,
    sigma_rv
)

end


"""
Build observation objects from mock StarData with noise injection.

For mock data fitting, this function adds Gaussian noise to the noiseless data
stored in StarData (generated by stardata_struct).

Parameters:
- star: Mock StarData object (noiseless observable data, stores uncertainties)
- epoch_mjd: Observation epoch in Modified Julian Date
- include_rv: Whether to build RV observation (default true for E and F)
- z_prior_sigma: If not nothing, build a PlanetZPriorObs with Normal(0, z_prior_sigma) in AU
- include_esc_vel: If true, build a PlanetEscVelObs using star.v2D / star.v2D_err [km/s]
- include_acc: If true, build a PlanetAccelObs using star.acc_ra / star.acc_dec [mas/yr²]

Returns:
- astrom: PlanetRelAstromObs (single-epoch relative position with added noise)
- pm: PlanetPMObs (single-epoch proper motion with added noise)
- acc: PlanetAccelObs (single-epoch acceleration with added noise)
- rv: PlanetRelativeRVObs or nothing (single-epoch RV with added noise, if available)
- zp: PlanetZPriorObs or nothing (LOS position prior)
- ev: PlanetEscVelObs or nothing (escape velocity constraint)
"""
function build_mock_observations(star::StarData, epoch_mjd::Float64;
                                  include_rv::Bool=true,
                                  z_prior_sigma::Union{Nothing,Float64}=nothing,
                                  include_esc_vel::Bool=false,
                                  include_acc::Bool=false,
                                )
   
    # 1. Single-epoch position with NOISE 
    astrom = PlanetRelAstromObs(
            (epoch = [epoch_mjd],
            ra    = [rand(Normal(star.ra,  star.sigma_ra))],
            dec   = [rand(Normal(star.dec, star.sigma_dec))],
            σ_ra  = [star.sigma_ra],
            σ_dec = [star.sigma_dec],
            cor   = [0.0]);
            name = "$(star.name)_pos"
        )
    
    # 2. Proper motion with NOISE 
    pm = PlanetPMObs(
        (epoch=[epoch_mjd],
         pmra=[rand(Normal(star.pm_ra, star.sigma_pm_ra))],
         pmdec=[rand(Normal(star.pm_dec, star.sigma_pm_dec))],
         σ_pmra=[star.sigma_pm_ra], σ_pmdec=[star.sigma_pm_dec], cor=[0.0]);
        name="$(star.name)_pm"
    )

    # 3. Acceleration with NOISE 
    acc = nothing
    if include_acc
        acc = PlanetAccelObs(
            (epoch = [epoch_mjd],
            accra  = [rand(Normal(star.acc_ra,  star.sigma_acc_ra))],
            accdec = [rand(Normal(star.acc_dec, star.sigma_acc_dec))],
            σ_accra  = [star.sigma_acc_ra],
            σ_accdec = [star.sigma_acc_dec],
            cor = [0.0]);
            name = "$(star.name)_acc"
        )
    end

    # 4. Radial velocity with NOISE 
    # Not substracting cluster RV here since stardata_struct already generates a peculiar RV relative to the cluster.
    rv = nothing
    if include_rv && !isnan(star.rv) && !isnan(star.rv_err)
        rv_noisy = rand(Normal(star.rv, star.rv_err))
        rv = PlanetRelativeRVObs(
        (epoch = epoch_mjd,
        rv    = rand(Normal(star.rv, star.rv_err)),
        σ_rv  = star.rv_err,
        variables = @variables begin end);
        name = "$(star.name)_rv"
        )
    end

    # 5. LOS position prior (z ~ Normal(0, σ_z) in AU)
    zp = nothing
    if z_prior_sigma !== nothing
        zp = PlanetZPriorObs(epoch_mjd, Normal(0.0, z_prior_sigma);
                              name="$(star.name)_zprior")
    end

    # 6. Escape velocity constraint
    ev = nothing
    if include_esc_vel
        ev = PlanetEscVelObs(epoch_mjd, star.v2D, star.v2D_err, v_esc_cluster_kms;
                              name="$(star.name)_escvel")
    end

    return astrom, pm, acc, rv, zp, ev
end
# ========== END MOCK MODIFICATIONS ==========

end

