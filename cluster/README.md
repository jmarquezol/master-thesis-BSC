# Wall-scan sweeps on the cluster

Each job walks a ladder of evolution times T and saves after every rung. If it hits the walltime,
resubmit the same script: it resumes from its checkpoint. To extend a finished ladder, raise the
Tmax at the bottom of its script and resubmit.

All jobs use the corrected transfer-matrix column (August 14 fix; the old column was missing a
memory channel at p != 0). Their results go to separate `*_bulk.jld2` caches. The superseded
scripts are in `archive/legacy_column/` — do not run them, and cancel any still in the queue.
The p = 0 job is unaffected by the fix.

These jobs continue the runs of 15-16 August rather than replacing them. The labels are unchanged,
so each ladder skips the rungs already in its cache and computes the half-integer ones between them.
Keep the checkpoints in `cluster/checkpoints/`; they are what makes the resume warm.

p = 1.5 has been dropped. It produced a single rung at T=2 in every route, one of them costing 42 h,
and its boundary dimension came out at 0.400 against 0.505 and 0.494 at p = 1.0. The scripts are in
`archive/dropped_p1.5/` if it is ever worth revisiting.

## Setup

```
git clone https://github.com/jmarquezol/master-thesis-BSC.git
cd master-thesis-BSC
module load julia/1.12.0
julia --project=. -e 'using Pkg; Pkg.instantiate()'
```

If you already have the repository: `git checkout main && git pull`. Do not update the packages —
`Manifest.toml` pins ITransverse to a specific commit.

## Submitting

From inside `cluster/`. The spectral arms carry the central charge at p != 0 and matter most if the
queue has to be cut short.

```
cd cluster
julia --project=.. wall_scan_cluster.jl preflight    # must print "preflight OK"

for p in 0.0 0.1 0.3 0.5 1.0; do
    sbatch submit_eigs_p${p}_bulk.slurm       # central charge from the eigenvalue phase
done

for p in 0.0 0.1 0.3 0.5 1.0; do
    sbatch submit_ent_p${p}_bulk.slurm        # temporal entropy
done

for p in 0.0 0.1 0.3 0.5 1.0; do
    sbatch submit_tower_p${p}_bulk.slurm      # boundary dimensions
done
```

## Notes

Three routes per coupling: `eigsweep` for the central charge from the eigenvalue phase, `entsweep`
for the temporal entropy, `towerscan` for the boundary dimensions.

`dT` must be a multiple of the Trotter step, since the chain holds T/dt + nbeta sites. dT=0.5 works
at dt=0.1, dT=0.25 needs dt=0.05. The driver refuses anything else.

Half-integer rungs matter more than reach. At p=0.3 an integer ladder gave c=0.41 with residuals
1.2e-3; the same range at dT=0.5 gave c=0.51 with residuals 3.7e-4.

Targets are Tmax 20 for the eigenvalue arms, 18 for the entropy arms and 8 for the towers, with p=0
at 24 on both routes. These are targets, not expectations: a walltime kill costs only the rung in
flight. The p=0 entropy arm runs longest because the finite-time correction converges slowly there
— fitting it over T<=14 gives c=0.578, over T<=18 gives 0.561, over T<=24 gives 0.514.

Which rungs are usable is decided in analysis, by the phase-advance criterion in
`analysis/session_scripts/eq3_fits.jl`, not by where a job happened to stop.

A checkpoint is only used to seed a longer chain than its own. Interleaving finer rungs below the
last completed one would otherwise abort the ladder, since a converged vector cannot be shrunk.

`WALL_RETRIES` (default 2) recomputes a rung from a fresh seed when its leading eigenvalue jumps
away from the previous one. Failures are seed-dependent, not a hard wall.

Check the thread scaling before filling 40 cores. Locally BLAS saturated at 2 threads and got slower
beyond it; if that holds here, concurrent jobs on fewer cores each will beat one job on 40.

## Where things land

- results: `results/data/cluster/` — this is the folder to send back
- checkpoints: `cluster/checkpoints/` — never delete these
- logs: `cluster/logs/`
