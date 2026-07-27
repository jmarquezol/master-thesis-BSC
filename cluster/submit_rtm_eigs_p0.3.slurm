#!/bin/bash
#SBATCH --job-name=wallscan_eigs_p0.3
#SBATCH --output=logs/wallscan_eigs_p0.3-%j.out
#SBATCH --error=logs/wallscan_eigs_p0.3-%j.err
#SBATCH --account=bsc21
#SBATCH --qos=gp_bsccase
#SBATCH --time=2-00:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=40

# EIGENVALUE-ONLY p-sweep arm: p=0.3, RTM route, Schur / eigvals-only mode, warm-started and
# checkpointed, T=2..14. Identical configuration to the earlier submit_rtm_p*.slurm full runs
# (chi=64, same strict cutoff schedule, same warm-start + checkpoint machinery) EXCEPT the block
# solver computes ONLY the leading spectrum and skips the eigenVECTOR work (temporal entropy +
# phase rigidity). The spectral observables we want from this — the dual-unitarity circle, the
# Eq.(3) central charge, the Eq.(4) boundary exponent, and the tower gaps — are all Rayleigh-
# quotient-like and survive the entanglement-barrier wall, so this arm is both cheaper per T and
# reaches larger T than the full-eigenvector runs (which stall once the eigenvector conditioning
# collapses). Output goes to its OWN cache, results/data/cluster/sweep_rtm_eigs_p0.3.jld2, so it
# never clobbers the full-eigenvector sweep_rtm_p0.3.jld2.
#
# If the job hits the walltime cap, just `sbatch` this same script again: it resumes from the last
# checkpointed T, truly warm (not just continuity-anchored).

cd "${SLURM_SUBMIT_DIR:-.}" || exit 1
mkdir -p logs

module load julia/1.12.0

export JULIA_NUM_THREADS=$SLURM_CPUS_PER_TASK
export OPENBLAS_NUM_THREADS=$SLURM_CPUS_PER_TASK

julia --project=. cluster/wall_scan_cluster.jl eigsweep 0.3 14
