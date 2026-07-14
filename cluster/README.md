# Master wall-scan sweep — BSC cluster quick-start

Three independent, full T=2..14 ladders of the same physics (Alcaraz ANNNI-type model, p=0.1,
quench to criticality), each isolating a different truncation knob:

| job | `trunc_mode` | notes |
|---|---|---|
| `submit_rtm.slurm` | `:rtm` (default) | re-derived χ=64 baseline, same config as the existing local run but extended to T=14 |
| `submit_rdm.slurm` | `:rdm` | truncates each vector independently — documented as better-conditioned near the degeneracy, but costlier |
| `submit_cutoff.slurm` | `:rtm`, tighter cutoff | `1e-10`/`1e-12` instead of the production `1e-8`/`1e-10` |

All three write into the **same** `results/data/nb13_wallscan_cluster.jld2` (keyed by
`(label, T)`), so they can run in any order, in parallel, and merge without conflict.

## Steps

1. **Clone and enter the repo:**
   ```
   git clone <repo-url>
   cd <repo>
   mkdir -p cluster/logs
   ```

2. **Instantiate the Julia environment — on a login node** (needs internet to fetch
   `ITransverse.jl` from GitHub; compute nodes are normally offline, so this must happen before
   any batch submission):
   ```
   julia --project=. -e 'using Pkg; Pkg.instantiate()'
   ```
   Sanity check it worked:
   ```
   julia --project=. -e 'using ITransverse, ITensorExpMPO; println("ok")'
   ```
   should print `ok` with no errors.

3. **Fill in the `<TODO>` placeholders** in all three `cluster/submit_*.slurm` files:
   - Which system are you on — **MN5** or **GRACE**? The `--qos` prefix differs (`gp_` on MN5,
     `ngp_` on GRACE).
   - `--account` — your BSC project code (e.g. `bsc21`).
   - `--qos` — the account's normal queue (e.g. `gp_bsccase`); there's usually also a `_debug`
     variant for short test jobs.
   - `--time` — a generous walltime. T=11..14 is uncharted for this project and likely the most
     expensive points yet. If a job hits the cap, it just stops — see step 6, nothing is lost.
   - `--cpus-per-task` — how many cores your account can request per node (MN5 general-purpose
     nodes have up to 112).
   - The `--output`/`--error`/`cd` absolute paths — point them at wherever you cloned the repo.
   - The `module load julia/...` line — check `module avail julia` for the exact name/version.

4. **Smoke-test one job first**, on the `_debug` qos with a short `--time` (e.g. `0:30:00`), before
   committing cluster allocation time to the full run:
   ```
   sbatch cluster/submit_rtm.slurm
   ```
   (with `--qos`/`--time` temporarily set to the debug/short values). Once it produces a sane
   `results/data/nb13_wallscan_cluster.jld2` for T=2, switch `--qos`/`--time` back to the real
   values for the actual submissions.

5. **Submit all three:**
   ```
   sbatch cluster/submit_rtm.slurm
   sbatch cluster/submit_rdm.slurm
   sbatch cluster/submit_cutoff.slurm
   ```
   They run independently — in parallel if the cluster has capacity.

6. **Monitor:**
   ```
   squeue -u $USER
   tail -f cluster/logs/wallscan_rtm-<jobid>.out      # etc. for rdm / cutoff
   ```
   Each per-T point logs a line like
   `[rtm64_full] T=9.0  stuck@1069  |θ0|=1.5443  peak=0.3331  r=[...]  4533s`
   as it completes.

7. **If a job hits its walltime cap**, it simply stops mid-ladder — no data is lost, since every
   completed T is saved to the cache as it finishes. Just resubmit the *same* `sbatch` command;
   it will skip every already-completed T and continue from where it left off.

8. **Send the results back** — once all three finish (or as far as they get):
   ```
   git add results/data/nb13_wallscan_cluster.jld2
   git commit -m "cluster master sweep results (rtm/rdm/cutoff, T=2..14)"
   git push
   ```
   If you'd rather hand the file back another way (scp, email), that's fine too — just let us know.

## What the output looks like

`results/data/nb13_wallscan_cluster.jld2` holds a `Dict` keyed by `(label, T)` — e.g.
`("rdm64", 9.0)` — where each value is a named tuple with `theta`, `theta_phys`, `s2_base` (the
trimmed Rényi-2 dome), `peak`, `rigidity`, `reason` (`"converged"` or `"stuck"`), `niters`, and
`elapsed`. Loadable with `JLD2.load(path, "done")`.
