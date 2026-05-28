#!/bin/bash
#SBATCH --account=def-vhenault
#SBATCH --job-name=mock_test_1
#SBATCH --nodes=1
#SBATCH --cpus-per-task=192
#SBATCH --mem-per-cpu=3G
#SBATCH --time=10:00:00
#SBATCH --output=/home/kenzhayd/projects/def-vhenault/kenzhayd/Ocen_IMBH_Mock_Analysis/configs/../mock_results/logs/mock_test_1_%j.out
#SBATCH --error=/home/kenzhayd/projects/def-vhenault/kenzhayd/Ocen_IMBH_Mock_Analysis/configs/../mock_results/logs/mock_test_1_%j.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=Mackenzie.hayduk@smu.ca

mkdir -p /home/kenzhayd/projects/def-vhenault/kenzhayd/Ocen_IMBH_Mock_Analysis/configs/../mock_results/logs
mkdir -p /home/kenzhayd/projects/def-vhenault/kenzhayd/Ocen_IMBH_Mock_Analysis/configs/../mock_results/run_outputs

module load julia/1.11.3

export JULIA_CONDAPKG_BACKEND=Null

julia --project=/home/kenzhayd/projects/def-vhenault/kenzhayd/octoIMBH_env -e 'using Pkg; Pkg.instantiate(); Pkg.add(["CairoMakie", "PairPlots", "Distributions", "Unitful", "UnitfulAstro", "Pigeons", "KernelDensity"])' 2>&1 | tee /home/kenzhayd/projects/def-vhenault/kenzhayd/Ocen_IMBH_Mock_Analysis/configs/../mock_results/logs/instantiate_${SLURM_JOB_ID}.log

julia --project=/home/kenzhayd/projects/def-vhenault/kenzhayd/octoIMBH_env -t 192 \
    /home/kenzhayd/projects/def-vhenault/kenzhayd/Ocen_IMBH_Mock_Analysis/launch_scripts/octo_orbit_direct_likelihoods_2.jl \
    /home/kenzhayd/projects/def-vhenault/kenzhayd/Ocen_IMBH_Mock_Analysis/configs/mock_default_2.toml 2>&1 | tee /home/kenzhayd/projects/def-vhenault/kenzhayd/Ocen_IMBH_Mock_Analysis/configs/../mock_results/logs/output_${SLURM_JOB_ID}.log
