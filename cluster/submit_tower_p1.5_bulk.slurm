#!/bin/bash
#SBATCH --job-name=tower_p1.5_bulk
#SBATCH --output=logs/tower_p1.5_bulk-%j.out
#SBATCH --error=logs/tower_p1.5_bulk-%j.err
#SBATCH --account=bsc21
#SBATCH --qos=gp_bsccase
#SBATCH --time=2-00:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=40

# Deep k=8 block at p=1.5; exploratory.
# SECOND TIER — submit only after the p <= 0.5 jobs are running. Corrected column
# (WALL_COLUMN=bulk5). The phase window shrinks as 1/v and v roughly doubles from p=0.5 to 1.5,
# the walltime, not Tmax, is the real cap: rungs stream to the cache as they finish, and a resubmission resumes warm from the checkpoint.
# WALL_DT below was set from the local Trotter check; do not change it independently of NBETA.

cd "${SLURM_SUBMIT_DIR:-.}" || exit 1
mkdir -p logs
module load julia/1.12.0
export JULIA_NUM_THREADS=$SLURM_CPUS_PER_TASK
export OPENBLAS_NUM_THREADS=$SLURM_CPUS_PER_TASK
export WALL_COLUMN=bulk5
export WALL_DT=0.025

julia --project=.. wall_scan_cluster.jl preflight || { echo "PREFLIGHT FAILED - bad environment, not starting the ladder"; exit 1; }

julia --project=.. wall_scan_cluster.jl towerscan 1.5 8
