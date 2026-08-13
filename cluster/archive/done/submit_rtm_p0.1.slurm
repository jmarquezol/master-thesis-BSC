#!/bin/bash
#SBATCH --job-name=rtm_p0.1
#SBATCH --output=logs/rtm_p0.1-%j.out
#SBATCH --error=logs/rtm_p0.1-%j.err
#SBATCH --account=bsc21
#SBATCH --qos=gp_bsccase
#SBATCH --time=2-00:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=40

# Full run (entropy + spectrum) at p=0.1, extending the ladder past T=8. The temporal entropy stops
# being meaningful beyond the entanglement wall at T~9, so this documents where the wall sits rather
# than producing usable entropies past it. Resubmit if it hits the walltime — it resumes warm.
# Writes to results/data/cluster/sweep_rtm_p0.1.jld2.

cd "${SLURM_SUBMIT_DIR:-.}" || exit 1
mkdir -p logs
module load julia/1.12.0
export JULIA_NUM_THREADS=$SLURM_CPUS_PER_TASK
export OPENBLAS_NUM_THREADS=$SLURM_CPUS_PER_TASK

julia --project=.. wall_scan_cluster.jl preflight || { echo "PREFLIGHT FAILED - bad environment, not starting the ladder"; exit 1; }

julia --project=.. wall_scan_cluster.jl psweep 0.1 12
