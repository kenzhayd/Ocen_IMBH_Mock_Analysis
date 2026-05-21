"""
    parse_config.jl — TOML configuration loader and prior parser

Provides helpers for reading a run configuration file and converting
prior specification strings into Distributions.jl / Octofitter objects.

Usage (from the main fitting script):
```julia
include(joinpath(@__DIR__, "parse_config.jl"))
cfg = load_config("configs/default.toml")
prior = parse_prior("Uniform(10, 2_000_000)")
```
"""

using TOML
using Distributions
using Octofitter: Sine, UniformCircular

# ── Config loading ──────────────────────────────────────────────────────

"""
    load_config(path::String) -> Dict

Load and return a TOML configuration file as a nested Dict.
"""
function load_config(path::String)
    isfile(path) || error("Configuration file not found: $path")
    return TOML.parsefile(path)
end

# ── Prior parsing ───────────────────────────────────────────────────────

# Helper: parse a number string that may contain underscores (e.g. "2_000_000")
_parse_num(s::AbstractString) = parse(Float64, replace(strip(s), "_" => ""))

"""
    parse_prior(spec::String) -> Distribution or Parameterization

Convert a prior specification string into a Distributions.jl object
or an Octofitter Parameterization (e.g. `UniformCircular`).

Supported forms:
- `"Uniform(lo, hi)"`
- `"Normal(mu, sigma)"`
- `"truncated(Normal(mu, sigma), lower=L)"`
- `"truncated(Normal(mu, sigma), lower=L, upper=U)"`
- `"Sine()"`
- `"UniformCircular()"`

Underscores in numbers are allowed (e.g. `"Uniform(10, 2_000_000)"`).
"""
function parse_prior(spec::String)
    s = strip(spec)

    # --- Parameterizations (not Distributions) ---
    s == "UniformCircular()" && return UniformCircular()

    # --- Simple zero-argument distributions ---
    s == "Sine()" && return Sine()

    # --- truncated(Normal(mu, sigma), lower=..., upper=...) ---
    m = match(r"^truncated\(\s*Normal\(\s*([^,]+),\s*([^)]+)\)\s*,(.+)\)$", s)
    if m !== nothing
        mu    = _parse_num(m[1])
        sigma = _parse_num(m[2])
        rest  = m[3]
        kw = Dict{Symbol, Float64}()
        ml = match(r"lower\s*=\s*([^\s,)]+)", rest)
        ml !== nothing && (kw[:lower] = _parse_num(ml[1]))
        mu_match = match(r"upper\s*=\s*([^\s,)]+)", rest)
        mu_match !== nothing && (kw[:upper] = _parse_num(mu_match[1]))
        return truncated(Normal(mu, sigma); kw...)
    end

    # --- Uniform(a, b) ---
    m = match(r"^Uniform\(\s*([^,]+),\s*([^)]+)\)$", s)
    if m !== nothing
        return Uniform(_parse_num(m[1]), _parse_num(m[2]))
    end

    # --- Normal(mu, sigma) ---
    m = match(r"^Normal\(\s*([^,]+),\s*([^)]+)\)$", s)
    if m !== nothing
        return Normal(_parse_num(m[1]), _parse_num(m[2]))
    end

    error("Unrecognized prior specification: \"$spec\". " *
          "Supported: Uniform, Normal, truncated(Normal(...), ...), Sine(), UniformCircular().")
end

# ── Companion prior lookup (with per-star overrides) ────────────────────

"""
    get_companion_prior(cfg, star_name::String, param::String) -> String

Return the prior specification string for `param` of star `star_name`.
Checks `cfg["priors"]["overrides"][star_name][param]` first; falls back
to `cfg["priors"]["companion_defaults"][param]`.
"""
function get_companion_prior(cfg, star_name::String, param::String)
    overrides = get(cfg, "priors", Dict())
    star_overrides = get(get(overrides, "overrides", Dict()), star_name, Dict())
    if haskey(star_overrides, param)
        return star_overrides[param]
    end
    return cfg["priors"]["companion_defaults"][param]
end

# ── Epoch helper ────────────────────────────────────────────────────────

# ── Data selection flags ───────────────────────────────────────────────

"""
    get_data_flag(cfg, star_name::String, data_type::String) -> Bool

Return whether observation type `data_type` should be included for star
`star_name`.  Checks `cfg["data"]["overrides"][star_name][data_type]`
first; falls back to `cfg["data"]["defaults"][data_type]`; defaults to
`true` if neither section exists.
"""
function get_data_flag(cfg, star_name::String, data_type::String)::Bool
    defaults  = get(get(cfg, "data", Dict()), "defaults", Dict())
    default_val = get(defaults, data_type, true)
    overrides = get(get(get(cfg, "data", Dict()), "overrides", Dict()), star_name, Dict())
    return Bool(get(overrides, data_type, default_val))
end

# ── Acceleration type selection ───────────────────────────────────────

"""
    get_accel_type(cfg, star_name::String) -> String

Return the acceleration likelihood type for `star_name`: `"vector"` (2D components,
default), `"magnitude"` (scalar |a_sky|), or `"none"` (excluded).

Checks `cfg["data"]["overrides"][star_name]["acceleration"]` first, then
`cfg["data"]["defaults"]["acceleration"]`; defaults to `"vector"` if absent.

For backwards compatibility `true` → `"vector"` and `false` → `"none"`.
"""
function get_accel_type(cfg, star_name::String)::String
    defaults    = get(get(cfg, "data", Dict()), "defaults", Dict())
    default_val = get(defaults, "acceleration", "vector")
    overrides   = get(get(get(cfg, "data", Dict()), "overrides", Dict()), star_name, Dict())
    val         = get(overrides, "acceleration", default_val)
    if val === true || val == "vector"
        return "vector"
    elseif val === false || val == "none" || val == false
        return "none"
    elseif val == "magnitude"
        return "magnitude"
    else
        error("Unrecognized acceleration type: $(repr(val)). " *
              "Expected \"vector\", \"magnitude\", true, or false.")
    end
end

# ── Z-prior configuration ─────────────────────────────────────────────

"""
    get_z_prior_sigma(cfg) -> Union{Nothing, Float64}

Return the σ_z value (in AU) for the LOS position prior, or `nothing`
if the `[data.z_prior]` section is absent.
"""
function get_z_prior_sigma(cfg)::Union{Nothing, Float64}
    z_cfg = get(get(cfg, "data", Dict()), "z_prior", nothing)
    z_cfg === nothing && return nothing
    return Float64(z_cfg["sigma_z_au"])
end

# ── Epoch helper ────────────────────────────────────────────────────────

"""
    get_epoch_mjd(cfg) -> Float64

Convert `cfg["epoch"]["year"]` to MJD via `Octofitter.years2mjd`.
"""
function get_epoch_mjd(cfg)
    return Octofitter.years2mjd(cfg["epoch"]["year"])
end
