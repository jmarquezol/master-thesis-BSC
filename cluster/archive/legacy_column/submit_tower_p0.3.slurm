#!/bin/bash
#SBATCH --job-name=tower_p0.3
#SBATCH --output=logs/tower_p0.3-%j.out
#SBATCH --error=logs/tower_p0.3-%j.err
#SBATCH --account=bsc21
#SBATCH --qos=gp_bsccase
#SBATCH --time=1-00:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=40

# k=8 eigenvalue ladder at p=0.3. The k-block test in notebook 8 only ever reached T=2, so the
# block-size question is open across the ladder rather than answered. towerscan already runs k=8
# with k_retry=10 and eigenvalues only, so no driver change is needed.
# Own cache: results/data/cluster/sweep_tower_p0.3.jld2.

cd "${SLURM_SUBMIT_DIR:-.}" || exit 1
mkdir -p logs
module load julia/1.12.0
export JULIA_NUM_THREADS=$SLURM_CPUS_PER_TASK
export OPENBLAS_NUM_THREADS=$SLURM_CPUS_PER_TASK

julia --project=.. wall_scan_cluster.jl preflight || { echo "PREFLIGHT FAILED - bad environment, not starting the ladder"; exit 1; }

julia --project=.. wall_scan_cluster.jl towerscan 0.3 10
