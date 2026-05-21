# Mock Data Fitting with Ocen_IMBH_analysis

## What Changed FROM Ocen_IMBH_analysis?

(octo_utils.jl):
- Added mock orbit generation: `make_star()`, `mock_data()`, `stardata_struct()`
- Added noise to mock data : `build_mock_observations()`

(octo_orbit_direct_likelihoods.jl):
- If `[mock]` section in config with `enabled=true`: uses mock data
- Otherwise: uses real data 

Initial config template:
- `configs/mock_default.toml` 

---

## How to Run Mock Fits Locally

### Default run locally 
```bash
julia --project=../../Octofitter_imbh.jl octo_orbit_direct_likelihoods.jl ../configs/mock_default.toml
```
### Default run on a interactive session
julia --project=../../octoIMBH_env octo_orbit_direct_likelihoods.jl ../configs/mock_default.toml


## File Organization

```
launch_scripts/
  ├── octo_utils.jl                          # Mock functions added 
  ├── octo_orbit_direct_likelihoods.jl       # Steps for mock fitting added 
  ├── plot_chain.jl                          # (unchanged)
  └── ...

configs/
  ├── default.toml                           # real data
  ├── mock_default.toml                      # mock data
  └── ...

Ocen_IMBH_analysis/
  ├── results/                               # Real fits
  └── mock_results/                          # Mock fits 



