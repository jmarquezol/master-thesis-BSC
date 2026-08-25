#!/bin/bash
#SBATCH --job-name=eigs_p0.5_a
#SBATCH --output=logs/eigs_p0.5_a-%j.out
#SBATCH --error=logs/eigs_p0.5_a-%j.err
#SBATCH --account=bsc21
#SBATCH --qos=gp_bsccase
#SBATCH --time=24:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=40

# Eigenvalue arm on the corrected column at p=0.5. The spectral route is the one that survives
# frustration: theta is a Rayleigh quotient, so it is accurate to second order in the eigenvector
# error, and it carries the central charge where the entropy route no longer can.
# dT<1 is deliberate. What limits the Eq.(3) fit is the number of points inside the window, not
# the reach in T: at p=0.3 an integer ladder gave c=0.41 with residuals 1.2e-3, and adding the
# half-integer rungs over the same range gave c=0.51 with residuals 3.7e-4.
# Own cache: data/cluster/sweep_rtm_eigs_p0.5_bulk.jld2

cd "${SLURM_SUBMIT_DIR:-.}" || exit 1
mkdir -p logs
module load julia/1.12.0
export JULIA_NUM_THREADS=$SLURM_CPUS_PER_TASK
export OPENBLAS_NUM_THREADS=$SLURM_CPUS_PER_TASK
export WALL_COLUMN=bulk5

julia --project=.. wall_scan_cluster.jl preflight || { echo "PREFLIGHT FAILED - bad environment, not starting the ladder"; exit 1; }

julia --project=.. wall_scan_cluster.jl eigsweep 0.5 20 1.0 6 rtm_eigs_p0.5_fine
