#!/bin/bash
#SBATCH --job-name=ksec_p0.1_plus
#SBATCH --output=logs/ksec_p0.1_plus-%j.out
#SBATCH --error=logs/ksec_p0.1_plus-%j.err
#SBATCH --account=bsc21
#SBATCH --qos=gp_bsccase
#SBATCH --time=2-00:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=40

# Eigenvalue ladder confined to the K=plus sector at p=0.1. The control arm: at p=0.1 the mixed fine
# ladder already settles at c = 0.42-0.44, so the sector ladders must reproduce it.
# inside a sector the displaced family does not exist, so the branch is the sector label and no
# pi-jump can occur. Rungs past the sector crossing may report "stuck"; the eigenvalue is still
# accurate there and the run resumes from its checkpoint on resubmission.
# Own cache: results/data/cluster/sweep_ksec_p0.1_plus.jld2.

cd "${SLURM_SUBMIT_DIR:-.}" || exit 1
mkdir -p logs
module load julia/1.12.0
export JULIA_NUM_THREADS=$SLURM_CPUS_PER_TASK
export OPENBLAS_NUM_THREADS=$SLURM_CPUS_PER_TASK

julia --project=.. wall_scan_cluster.jl preflight || { echo "PREFLIGHT FAILED - bad environment, not starting the ladder"; exit 1; }

julia --project=.. wall_scan_cluster.jl ksector 0.1 plus 18
