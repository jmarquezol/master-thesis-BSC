#!/bin/bash
#SBATCH --job-name=tower_p0.0_bulk
#SBATCH --output=logs/tower_p0.0_bulk-%j.out
#SBATCH --error=logs/tower_p0.0_bulk-%j.err
#SBATCH --account=bsc21
#SBATCH --qos=gp_bsccase
#SBATCH --time=24:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=40

# Reconstructed. This arm ran in an earlier round and its submission file was not kept; the
# script below is the p=0.3 one with the coupling changed, and its arguments reproduce the
# shipped cache (T=2 to 8 in steps of 0.5, complete).
#
# Deep k=8 block on the corrected column at p=0.0, for the tower figure.
# The chain is nearest-neighbour only here, so the three-site column would already be
# exact, but this arm was run on the five-site one like the others.
# WALL_COLUMN=bulk5 builds the evolution operator on five sites and takes the true middle
# tensor, giving temporal site dimension 13 instead of 7, at roughly 2-4x the cost per
# iteration. The driver appends _bulk to the label, so this can never overwrite a legacy cache.
# Own cache: data/cluster/sweep_tower_p0.0_bulk.jld2

cd "${SLURM_SUBMIT_DIR:-.}" || exit 1
# the scripts live one level down, so step back up if sbatch was run from there
[ -f wall_scan_cluster.jl ] || cd .. || exit 1
mkdir -p logs
module load julia/1.12.0
export JULIA_NUM_THREADS=$SLURM_CPUS_PER_TASK
export OPENBLAS_NUM_THREADS=$SLURM_CPUS_PER_TASK
export WALL_COLUMN=bulk5

julia --project=.. wall_scan_cluster.jl preflight || { echo "PREFLIGHT FAILED - bad environment, not starting the ladder"; exit 1; }

julia --project=.. wall_scan_cluster.jl towerscan 0.0 8 0.5
