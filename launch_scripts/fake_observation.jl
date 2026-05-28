using TOML
using Distributions
using Octofitter
using OctofitterRadialVelocity

_obs_scalar(obs, col::Symbol) = getproperty(obs.table, col)[1]

function save_observations_toml(path, star_names, epoch_mjd,
                                astrom_obs, pm_obs, acc_obs, rv_obs)
    data = Dict{String, Any}(
        "epoch_mjd" => epoch_mjd,
        "stars" => Dict{String, Any}(),
    )

    for name in star_names
        star_data = Dict{String, Any}(
            "ra" => _obs_scalar(astrom_obs[name], :ra),
            "dec" => _obs_scalar(astrom_obs[name], :dec),
            "sigma_ra" => _obs_scalar(astrom_obs[name], Symbol("σ_ra")),
            "sigma_dec" => _obs_scalar(astrom_obs[name], Symbol("σ_dec")),

            "pmra" => _obs_scalar(pm_obs[name], :pmra),
            "pmdec" => _obs_scalar(pm_obs[name], :pmdec),
            "sigma_pmra" => _obs_scalar(pm_obs[name], Symbol("σ_pmra")),
            "sigma_pmdec" => _obs_scalar(pm_obs[name], Symbol("σ_pmdec")),

            "has_acceleration" => acc_obs[name] !== nothing,
            "has_rv" => rv_obs[name] !== nothing,
        )

        if acc_obs[name] !== nothing
            star_data["accra"] = _obs_scalar(acc_obs[name], :accra)
            star_data["accdec"] = _obs_scalar(acc_obs[name], :accdec)
            star_data["sigma_accra"] = _obs_scalar(acc_obs[name], Symbol("σ_accra"))
            star_data["sigma_accdec"] = _obs_scalar(acc_obs[name], Symbol("σ_accdec"))
        end

        if rv_obs[name] !== nothing
            star_data["rv"] = _obs_scalar(rv_obs[name], :rv)
            star_data["sigma_rv"] = _obs_scalar(rv_obs[name], Symbol("σ_rv"))
        end

        data["stars"][name] = star_data
    end

    open(path, "w") do io
        TOML.print(io, data)
    end
end

function load_observations_toml(path, star_names, epoch_mjd;
                                z_prior_sigma::Union{Nothing, Float64}=nothing)
    data = TOML.parsefile(path)

    astrom_obs = Dict{String, Any}()
    pm_obs     = Dict{String, Any}()
    acc_obs    = Dict{String, Any}()
    rv_obs     = Dict{String, Any}()
    zp_obs     = Dict{String, Any}()
    ev_obs     = Dict{String, Any}()

    for name in star_names
        d = data["stars"][name]

        astrom_obs[name] = PlanetRelAstromObs(
            (epoch=[epoch_mjd],
             ra=[Float64(d["ra"])],
             dec=[Float64(d["dec"])],
             σ_ra=[Float64(d["sigma_ra"])],
             σ_dec=[Float64(d["sigma_dec"])],
             cor=[0.0]);
            name="$(name)_pos"
        )

        pm_obs[name] = PlanetPMObs(
            (epoch=[epoch_mjd],
             pmra=[Float64(d["pmra"])],
             pmdec=[Float64(d["pmdec"])],
             σ_pmra=[Float64(d["sigma_pmra"])],
             σ_pmdec=[Float64(d["sigma_pmdec"])],
             cor=[0.0]);
            name="$(name)_pm"
        )

        acc_obs[name] = if get(d, "has_acceleration", false)
            PlanetAccelObs(
                (epoch=[epoch_mjd],
                 accra=[Float64(d["accra"])],
                 accdec=[Float64(d["accdec"])],
                 σ_accra=[Float64(d["sigma_accra"])],
                 σ_accdec=[Float64(d["sigma_accdec"])],
                 cor=[0.0]);
                name="$(name)_acc"
            )
        else
            nothing
        end

        rv_obs[name] = if get(d, "has_rv", false)
            PlanetRelativeRVObs(
                (epoch=epoch_mjd,
                 rv=Float64(d["rv"]),
                 σ_rv=Float64(d["sigma_rv"]));
                name="$(name)_rv",
                variables=@variables begin end
            )
        else
            nothing
        end

        zp_obs[name] = z_prior_sigma === nothing ? nothing :
            PlanetZPriorObs(epoch_mjd, Normal(0.0, z_prior_sigma);
                             name="$(name)_zprior")

        ev_obs[name] = nothing
    end

    return astrom_obs, pm_obs, acc_obs, rv_obs, zp_obs, ev_obs
end