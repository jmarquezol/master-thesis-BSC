# Wall-scan sweeps on the cluster

Each job walks a ladder of evolution times T and saves after every rung. If it hits the walltime,
resubmit the same script: it resumes from its checkpoint. To extend a finished ladder, raise the
Tmax at the bottom of its script and resubmit.

All jobs use the corrected transfer-matrix column (August 14 fix; the old column was missing a
memory channel at p != 0). Their results go to separate `*_bulk.jld2` caches. The superseded
scripts are in `archive/legacy_column/` — do not run them, and cancel any still in the queue.
The p = 0 job is unaffected by the fix.

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

From inside `cluster/`, in this order. The queue was reorganised on 2026-08-16 after the local runs
showed which route survives frustration; see "What changed" below before resubmitting anything.

```
cd cluster
julia --project=.. wall_scan_cluster.jl preflight    # must print "preflight OK"

sbatch submit_eigs_p0.1_bulk.slurm        # spectral arms: the route that works at p != 0
sbatch submit_eigs_p0.3_bulk.slurm
sbatch submit_eigs_p0.5_bulk.slurm
sbatch submit_eigs_p1.0_bulk.slurm
sbatch submit_eigs_p1.5_bulk.slurm

sbatch submit_ent_p0.0_bulk.slurm         # entropy arms: p=0 carries the absolute result
sbatch submit_ent_p0.1_bulk.slurm
sbatch submit_ent_p0.3_bulk.slurm
sbatch submit_ent_p0.5_bulk.slurm

sbatch submit_rtm_eigs_p0.0_fine.slurm    # p=0 spectral control

sbatch submit_tower_p0.1_bulk.slurm       # k=8 blocks for the boundary tower
sbatch submit_tower_p0.3_bulk.slurm
sbatch submit_tower_p0.5_bulk.slurm
sbatch submit_tower_p1.0_bulk.slurm
sbatch submit_tower_p1.5_bulk.slurm
```

## What changed on 2026-08-16, and why

**Added: five spectral arms at p != 0.** The queue had none. The eigenvalues are Rayleigh quotients,
accurate to second order in the eigenvector error, so the spectral route survives frustration where
the entropy route does not. It is what measures the central charge at p != 0.

**Every sweep now runs at dT = 0.5 or finer.** What limits these fits is the number of points inside
the window, not the reach in T. At p=0.3 an integer ladder gave c = 0.41 with residuals 1.2e-3;
adding half-integer rungs over the same range gave c = 0.51 with residuals 3.7e-4. Densify, do not
extend.

**Tmax cut to just past the usable window.** Measured locally, the last trustworthy rung is T ~ 20
at p=0, 11 at p=0.1, 7 at p=0.3, 4 at p=0.5. Rungs past that come back broken, so the walltime they
used is better spent inside the window.

**The p = 1.0 and 1.5 jobs keep headroom above that**, because the conformal predictions are
asymptotic in T and those windows are short in absolute terms. The data say short windows there are
not a problem in themselves: at T=2 the measured x1 misses 1/2 by 2.3% at p=0.1, 0.5% at p=0.3 and
0.3% at p=0.5, so at fixed T the more frustrated couplings are the more conformal ones, which is
what the rising velocity would predict. But the one rung we have at p=1.0 sits at T=1.2 and is the
worst of the set, so there is a floor in absolute T. Do not shorten these two further.

**Superseded (moved to `archive/superseded_2026-08-16/`).** The entropy arms at p = 1.0 and 1.5:
their window closes around T ~ 2-3, so the arm returns almost nothing usable. The four `_dt0.1`
Trotter twins: the norm drift after eight VD2 layers is already 1.56 and 6.96 at those couplings, so
they would spend two days confirming that the standard step is inadequate there.

**Driver: `towerscan` now takes a dT argument**, and every sweep mode checks it. Only multiples of
the Trotter step are realisable, because the chain holds T/dt + nbeta sites: a rung asked for T=2.25
at dt=0.1 gets 22 sites and actually runs at 2.2. Mixing two T grids scatters every phase-derived
read. **dT = 0.5 is safe at dt = 0.1; dT = 0.25 needs dt = 0.05.** The driver now refuses the bad
combinations rather than running them silently.

## Where things land

- results: `results/data/cluster/` — this is the folder to send back
- checkpoints: `cluster/checkpoints/` — never delete these
- logs: `cluster/logs/`
