#!/bin/bash
#SBATCH --chdir .
#SBATCH --job-name=wallscan_rdm-%j
#SBATCH --output=<TODO: absolute path to your cloned repo>/cluster/logs/wallscan_rdm-%j.out
#SBATCH --error=<TODO: absolute path to your cloned repo>/cluster/logs/wallscan_rdm-%j.err
#SBATCH --account=<TODO: your BSC project code, e.g. bsc21>
#SBATCH --qos=<TODO: e.g. gp_bsccase on MN5, ngp_bsccase on GRACE — confirm which system first;
#              use the _debug variant (e.g. gp_debug) for a quick smoke test before the real run>
#SBATCH --time=<TODO: DD-HH:MM:SS, e.g. 2-00:00:00 — :rdm truncates each vector independently and
#              is documented as costlier than :rtm; request generously, the job is resumable>
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=<TODO: e.g. 112 on an MN5 general-purpose node — check what your account
#              can request>

# The RDM truncation mode (trunc_mode=:rdm), full T=2..14 ladder — tests whether the wall (found
# gap-limited at T=9->10 for :rtm at both χ=32 and χ=64) moves when the truncation stays
# well-conditioned through the near-degeneracy instead of relying on the (potentially ill-
# conditioned there) joint RTM SVD. Run alongside submit_rtm.slurm and submit_cutoff.slurm.

module purge
module load <TODO: BSC's julia module + version, e.g. `module load julia/1.12` — check
#           `module avail julia`>

export JULIA_NUM_THREADS=$SLURM_CPUS_PER_TASK
export OPENBLAS_NUM_THREADS=$SLURM_CPUS_PER_TASK

cd <TODO: absolute path to your cloned repo>
julia --project=. cluster/wall_scan_cluster.jl rdm
