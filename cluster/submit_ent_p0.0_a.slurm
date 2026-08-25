#!/bin/bash
#SBATCH --job-name=ent_p0.0_a
#SBATCH --output=logs/ent_p0.0_a-%j.out
#SBATCH --error=logs/ent_p0.0_a-%j.err
#SBATCH --account=bsc21
#SBATCH --qos=gp_bsccase
#SBATCH --time=24:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=40

# Entropy arm on the corrected column at p=0.3 — the coupling where the legacy dome never read conformal.
# CORRECTED COLUMN (WALL_COLUMN=bulk5): temporal site dimension 13 instead of 7, roughly 2-4x the
# cost per iteration of the legacy runs. Tmax is deliberately modest for this first calibration
# round; resubmit with a larger Tmax to extend — the ladder resumes from its own checkpoint.
# Own cache: results/data/cluster/sweep_ent_p0.0_bulk.jld2 — never shared with a legacy cache.

cd "${SLURM_SUBMIT_DIR:-.}" || exit 1
mkdir -p logs
module load julia/1.12.0
export JULIA_NUM_THREADS=$SLURM_CPUS_PER_TASK
export OPENBLAS_NUM_THREADS=$SLURM_CPUS_PER_TASK
export WALL_COLUMN=bulk5

julia --project=.. wall_scan_cluster.jl preflight || { echo "PREFLIGHT FAILED - bad environment, not starting the ladder"; exit 1; }

julia --project=.. wall_scan_cluster.jl entsweep 0.0 24 1.0 11.5 ent_p0.0
