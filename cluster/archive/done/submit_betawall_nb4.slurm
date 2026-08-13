#!/bin/bash
#SBATCH --job-name=betawall_nb4
#SBATCH --output=logs/betawall_nb4-%j.out
#SBATCH --error=logs/betawall_nb4-%j.err
#SBATCH --account=bsc21
#SBATCH --qos=gp_bsccase
#SBATCH --time=2-00:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16

# Control arm for the beta0-vs-wall test: p=0.1 at nbeta=4 (beta0=0.2, our default), T up to 14.
# Same code path as the nbeta=12 arm, so the comparison has no configuration confound. It also
# extends the existing warm p=0.1 baseline, which currently stops at T=8 -- right at the wall.

cd "$(dirname "$0")" || exit 1
mkdir -p logs
module load julia/1.12.0
export JULIA_NUM_THREADS=$SLURM_CPUS_PER_TASK
export OPENBLAS_NUM_THREADS=$SLURM_CPUS_PER_TASK
julia --project=.. wall_scan_cluster.jl preflight || { echo "PREFLIGHT FAILED - bad environment, not starting the ladder"; exit 1; }

julia --project=.. wall_scan_cluster.jl betawall 0.1 4 14
