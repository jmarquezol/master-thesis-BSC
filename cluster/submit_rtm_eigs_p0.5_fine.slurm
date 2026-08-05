#!/bin/bash
#SBATCH --job-name=eigs_p0.5_fine
#SBATCH --output=logs/eigs_p0.5_fine-%j.out
#SBATCH --error=logs/eigs_p0.5_fine-%j.err
#SBATCH --account=bsc21
#SBATCH --qos=gp_bsccase
#SBATCH --time=2-00:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=40

# half-integer rungs at p=0.5 (T=2.5,3.5,...) to resolve the phase winding.
# The rungs merge into the same cache as the unit ladder, so the analysis sees one denser ladder.
# Resubmit this script if it hits the walltime — it resumes from the last completed rung.
# Writes to results/data/cluster/sweep_rtm_eigs_p0.5.jld2.

cd "${SLURM_SUBMIT_DIR:-.}" || exit 1
mkdir -p logs
module load julia/1.12.0
export JULIA_NUM_THREADS=$SLURM_CPUS_PER_TASK
export OPENBLAS_NUM_THREADS=$SLURM_CPUS_PER_TASK

julia --project=.. wall_scan_cluster.jl preflight || { echo "PREFLIGHT FAILED - bad environment, not starting the ladder"; exit 1; }

julia --project=.. wall_scan_cluster.jl eigsweep 0.5 20 0.5
