#!/bin/bash
#SBATCH --job-name=eigs_p0.0_a
#SBATCH --output=logs/eigs_p0.0_a-%j.out
#SBATCH --error=logs/eigs_p0.0_a-%j.err
#SBATCH --account=bsc21
#SBATCH --qos=gp_bsccase
#SBATCH --time=2-00:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=40

# Half-integer rungs at the integrable point, the control the fine ladders never had.
# The fine and unit ladders disagree by 20% on c at p=0.1 and only the fine one settles, so the
# same analysis has to be run where the answer is known. Resubmit if it hits the walltime — it
# resumes from the last completed rung.
# Own cache: data/cluster/sweep_eigs_p0.0_bulk.jld2.

cd "${SLURM_SUBMIT_DIR:-.}" || exit 1
mkdir -p logs
module load julia/1.12.0
export JULIA_NUM_THREADS=$SLURM_CPUS_PER_TASK
export OPENBLAS_NUM_THREADS=$SLURM_CPUS_PER_TASK

julia --project=.. wall_scan_cluster.jl preflight || { echo "PREFLIGHT FAILED - bad environment, not starting the ladder"; exit 1; }

julia --project=.. wall_scan_cluster.jl eigsweep 0.0 24 1.0 12.5 rtm_eigs_p0.0_fine
