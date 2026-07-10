# Open question: the ± eigenvalue pairing in the Alcaraz transfer-matrix spectrum

Status: OPEN. Not part of the closed barrier campaign (§18 of CLAUDE.md) — this is a new thread
that surfaced while explaining notebook 5's "the gap closes" section, and it is not derived or
explained anywhere in the current codebase or CLAUDE.md. This document records what is confirmed
numerically, candidate physical explanations (all speculative, ranked by plausibility), and concrete
next steps to disambiguate them.

---

## 1. The empirical finding

Reading the leading `k=4` transfer-matrix eigenvalues from the converged master sweep
(`results/data/nb8_master.jld2`, produced by NB7's block-PM ladder) shows a robust pattern for the
frustrated Alcaraz model (`p=0.1`):

- At every `T` we checked (2, 3, 4, 5, 6, 7, 8), there is a "leading pair" — the physical `λ0`
  (`i0`, selected by `pick_phys`'s complex-value continuity) and a partner eigenvalue whose phase is
  ≈π away from `λ0`'s and whose modulus is close to `λ0`'s.
- This pairing is **not exact** and **not static**: at `T=2` the modulus ratio is only 0.94 (14%
  off from a perfect pair) and the phase offset is −2.87 rad (0.27 rad off from exactly π). By
  `T=8` the modulus ratio has tightened to 0.995–1.004 and the phase offset to within ~0.07 rad of
  π. The pairing **sharpens monotonically as T grows** — it is not a fixed algebraic symmetry
  present from the start, it is an emergent, asymptotically-tightening near-degeneracy.
- At small `T` (`T=2`) only the single leading pair is clean; the other two `k=4` Ritz values are
  unpaired (modulus ratios 0.65–0.89, phases not near π apart). By `T=5`–`8`, more of the `k=4`
  eigenvalues satisfy the pairing criterion simultaneously — the spectrum organizes into two
  sub-clusters ≈π apart. In other words: **the leading pair is the first instance of the same wall
  phenomenon** already documented in CLAUDE.md §17/§18 (the widening near-degenerate band that
  destroys eigenvector conditioning) — it is not a separate effect from the barrier, it looks like
  its earliest signature.

**The critical control experiment: Ising (`p=0`) shows none of this.** Running the identical
pairwise check on the `p=0` sweep — same cache, same asymmetric VD2 exponentiation kernel, same
`block_transfer_eigs` oblique-Rayleigh-Ritz solver, only the coupling coefficients differ (`p=0`
simply zeroes the NNN and `XX` entries of the same 5-state FSM described in CLAUDE.md §3) — finds
**zero pairs at any T=2, 5, 8**. Ising's four leading phases stay clustered together; none sit ≈π
apart at any point in this range.

This is the crux fact any explanation has to account for: **the pairing is not an artifact of the
transverse/VD2 construction** (both models are built by the exact same machinery), so it must come
from physics that is present at `p≠0` and absent at `p=0`.

### What is already ruled out

- **Not the global spin-flip symmetry sector structure.** `H` commutes with `P=∏_i σ^x_i` for
  *every* `p` (all four Hamiltonian terms — `ZZ`, `p·ZZ`(NNN), `λX`, `pλ·XX` — are invariant under
  simultaneously flipping every `σ^z_i → −σ^z_i`), so `P` cannot be what distinguishes `p=0` from
  `p≠0`. Moreover CLAUDE.md §10/§17 already established that `|X+⟩` is purely `P`-even and that the
  Alcaraz near-degeneracy lives *within* that even sector — sector projection was already tried and
  explicitly found not to help. So this is a different, subtler within-sector phenomenon, not the
  same mechanism as XXZ's exact Néel-quench `P`-sector degeneracy (§18), even though both end up
  looking like "an extra near-degenerate partner."
- **Not generic realness of H.** `H` being a real symmetric matrix guarantees `U(dt)=\exp(-iHdt)`'s
  eigenvalues come in complex-conjugate pairs `(λ, λ^*)` — but that is a *different* pairing
  (reflection across the real axis) than the observed `(λ, −λ)`-type pairing (reflection through the
  origin, ≈π phase separation). Realness alone does not produce the latter.

---

## 2. Candidate physical mechanisms (speculative — none confirmed)

Ranked roughly by how well each survives the `p=0` control test.

### (a) Loss of integrability under Jordan–Wigner (most promising, least developed)
At `p=0` the model is exactly the standard transverse-field Ising chain (CLAUDE.md §2: "`p=0`
exactly recovers the standard integrable TFIM"), solvable by Jordan–Wigner as free fermions. Any
coupling beyond nearest-neighbour — in particular the `p·ZZ` NNN term — generically becomes a
**quartic (interacting) fermion term** under the same transformation, breaking free-fermion
integrability for any `p≠0`. This is consistent with the pairing turning on exactly where
integrability breaks. It is a plausible *correlation* but not yet a *mechanism*: it is not obvious
a priori why "the fermion theory becomes interacting" should produce a transfer-matrix spectrum with
eigenvalues pairing at ≈π phase separation specifically. This would need the explicit JW-transformed
Hamiltonian for general `p` worked out, and a check of whether analogous eigenvalue-pairing shows up
in known results for interacting/frustrated Ising-type Floquet or transfer-matrix spectra in the
literature.

### (b) An unidentified residual symmetry of the asymptotic (dual-unitary) transfer operator
CLAUDE.md's own findings (NB6, NB9, §11 Step 6) establish that the transfer matrix approaches a
*rescaled unitary* operator as `T→∞` (emergent dual unitarity — `λ0(T)` traced in the complex plane
lies on a circle of near-constant radius). The pairing's key qualitative feature — it **sharpens
monotonically with T** rather than being exact from `T=0` — is exactly the signature you'd expect if
it is tied to this same large-`T` approach: if the asymptotic unitary transfer operator possesses
some additional discrete (Z2-like) symmetry, nearby-but-not-yet-converged eigenvalues could
increasingly organize into a symmetric pattern as the operator gets closer to its unitary fixed
point.

**Tension with the control experiment**: Ising is the *cleanest* realization of dual unitarity in
this project (CLAUDE.md §17: symmetric Murg + Takagi reaches `T=14` with a long clean window,
markedly better-conditioned than Alcaraz's asymmetric route) — if approach-to-dual-unitarity were
the generic driver of ± pairing, Ising should show it *most* clearly, not not at all. A possible
resolution: Ising's `p=0` run above uses the *same asymmetric* VD2 + oblique block-PM pipeline as
Alcaraz (not NB6's symmetric Murg + Takagi route), so this is a fair apples-to-apples comparison —
the absence of pairing is not a solver artifact. If this hypothesis is right, the residual symmetry
would have to be one that Ising's asymptotic unitary transfer operator *lacks* and Alcaraz's *has* —
plausibly tied back to hypothesis (a): the extra symmetry may only be present/relevant once the
fermion theory is interacting.

### (c) Self-duality of H at λ=1
`H(λ,p) = λ·H(1/λ,p)` holds for *every* `p` (CLAUDE.md §2), so self-duality alone can't distinguish
`p=0` from `p≠0` either — at the critical point `λ=1` it's a trivial statement (`H(1,p)=H(1,p)`) for
both models. Included for completeness but currently the weakest candidate; would only become
relevant if the pairing's *strength* (not just presence/absence) tracked distance from `λ=1` in a
way that couples to `p` — untested (see Next Steps §3.4).

### (d) A bipartite/chiral sublattice structure from mixing NN and NNN couplings
On a chain, NN bonds `(i,i+1)` always connect opposite-parity sites (even–odd), while NNN bonds
`(i,i+2)` always connect same-parity sites (even–even or odd–odd). Two-coloring the chain by site
parity, the pure-NN Ising Hamiltonian (`p=0`, plus the onsite field) has a mix of inter-sublattice
(`ZZ`) and onsite (`λX`) terms; turning on `p` *adds* a purely intra-sublattice term (`p·ZZ`,NNN)
plus another inter-sublattice term (`pλ·XX`). A matrix with *purely* off-diagonal block structure
between two equal-size subspaces provably has an exactly `λ↔−λ` symmetric spectrum (the standard
chiral/bipartite-lattice mechanism, e.g. particle–hole symmetry in tight-binding models). This is
an attractive, easily-stated mechanism, but **it does not obviously survive contact with the
model**: the onsite field term `λX_i` is diagonal in the sublattice decomposition (a "same-site"
term, neither purely inter- nor intra-), which breaks the clean off-diagonal block structure the
theorem needs — for *either* `p=0` or `p≠0`. It is not clear why adding the NNN/XX terms would
create the missing chiral structure rather than further breaking it. Would need to be checked
directly on a small explicit matrix rather than argued qualitatively (see Next Steps §3.5).

### (e) Boundary-state visibility, not spectrum existence (reframes the question)
An alternative framing: maybe Ising's *full* transfer matrix spectrum *also* has a `−λ0`-type
partner somewhere, and the true difference is not "Ising's spectrum lacks the pairing" but "`|X+⟩`
has negligible overlap with Ising's partner eigenvector, so it never appears among the leading `k=4`
Ritz values reached by the block power method," whereas Alcaraz's `|X+⟩` boundary happens to couple
comparably to both members of the pair. This would mean the physically interesting object is not
"why does `p≠0` create new eigenvalues" but "why does the same boundary state excite a pre-existing
pair unevenly depending on `p`" — a genuinely different (and more tractable) question. Directly
testable by increasing `k` in the Ising block PM (see Next Steps §3.2).

---

## 3. Concrete next steps to disambiguate

Ordered roughly by expected information gained per unit of compute — cheap diagnostics first.

### 3.1 Fine `p`-sweep near zero
Run the same pairwise-eigenvalue check (already scripted, see the transcript of this
investigation) at `p = 0.01, 0.02, 0.05, 0.1` (fixed `T`, say `T=5`). Does the pairing strength
(modulus-ratio deviation from 1, phase-offset deviation from π) turn on **continuously** from `p=0`
(consistent with a smooth integrability-breaking mechanism, hypothesis (a)) or does it appear
**abruptly** past some threshold `p*` (consistent with a bifurcation / level-crossing mechanism)?
Cheap: reuses the existing `nb8_master.jld2`-style block-PM machinery, no new algorithm needed.

### 3.2 Increase `k` for Ising specifically
Re-run `block_transfer_eigs` for `p=0` with `k=8` or `k=16` instead of `4`. If a `−λ0`-like
eigenvalue eventually appears further down the spectrum (just with far smaller `k=4`-invisible
weight), that supports hypothesis (e) — the pairing is a property of the full transfer matrix
regardless of `p`, only its *visibility* to `|X+⟩` differs. If no such eigenvalue appears even at
`k=16`, that is stronger evidence the pairing is genuinely `p≠0`-specific physics (hypotheses
(a)/(b)/(d)).

### 3.3 Try a different (symmetry-breaking) boundary state for Ising
Re-run the `p=0` sweep with the fixed-BC boundary state `|Z+⟩` instead of `|X+⟩` (already a
supported `init_state` per CLAUDE.md §13 — `|Z±⟩` is the "fixed BC, slow convergence,
`x₁=2`" alternative quench used as a cross-check elsewhere in the project). If pairing appears for
Ising under this different boundary condition, that strongly supports the boundary-state-visibility
reframing (e) over a `p≠0`-specific mechanism.

### 3.4 λ-sweep away from criticality, fixed p
Fix `p=0.1` and sweep `λ` away from `1` (e.g. `λ = 0.8, 1.0, 1.2`) at fixed `T`. Does the pairing
weaken away from the self-dual/critical point (supporting hypothesis (c), or more generally tying
the pairing to criticality itself) or is it insensitive to `λ` (ruling out both self-duality and a
purely-critical-point origin)?

### 3.5 Explicit small-system chiral-structure check
For a short chain (small enough for exact diagonalization or a small-bond-dimension exact transfer
matrix — reuse NB11's exact dense-diagonalization ground truth, already built and validated in the
barrier campaign), attempt to construct an explicit sublattice/parity operator `S` (diagonal,
`±1`-valued on some natural basis of the rotated temporal physical index) and check numerically
whether `S T S^{-1} = -T` holds exactly, approximately, or not at all, for `p=0` vs `p≠0`. This is
the most rigorous, "prove it algebraically" version of hypothesis (d), and NB11's existing exact
ground truth infrastructure means this doesn't require new solver code — only a new symmetry check
on matrices that are already being computed there.

### 3.6 Cross-check against the dual-unitarity circle — DONE (2026-07-09), result: NOT supported

**Test performed.** Using only the already-cached `nb8_master.jld2` (`p=0.1`, no new computation),
compared two things across the pre-wall window `T=3..9`: (i) the pairing tightness — `dev_mod =
|modulus_ratio-1|` and `dev_phase = |phase_offset-π|` for the `(i0, partner)` pair, partner found by
the same pairwise-search criterion used in §1; (ii) a convergence proxy for `λ0(T)` itself —
step-to-step motion `|λ0(T)-λ0(T-1)|`, chosen over "deviation from an assumed asymptote" because the
raw `|λ0(T)|` series is *not* monotonically settling within `T≤12` (it decreases `T=2→9` then rises
again `T=9→12`, the latter already past the wall), so any fixed reference point would be arbitrary
and the `T≥10` points are contaminated anyway.

**Result — the two do NOT track together.** `dev_phase` drops sharply and monotonically over
`T=3..9` (ratio to its `T=3` value: `1.00 → 0.75 → 0.60 → 0.51 → 0.43 → 0.38 → 0.34`) — the pairing
is genuinely tightening. Over the *same* window, `|λ0(T)-λ0(T-1)|` stays essentially flat (`1.99 →
1.99 → 1.98 → 1.98 → 1.98 → 1.98 → 1.97`, ratio to `T=3` never leaving `[0.99,1.00]`) — `λ0` is
sweeping its circle at essentially constant angular speed, showing no sign of decelerating/settling
over this range at all. (At `T=10`, both quantities jump — `dev_phase` back up to `0.91`×ref,
`|Δλ0|` up to `1.12`×ref — consistent with `T=10` being past the wall, excluded from the comparison.)

**Reading.** This is a clean negative result for hypothesis (b) as stated: if the pairing were a
signature of `λ0` itself settling onto its asymptotic circle, the two should decay on comparable
timescales — they visibly don't. `λ0`'s motion is already essentially steady from `T=3` (consistent
with `λ0` having *already* found its circular orbit early), while the pairing keeps sharpening
*independently* over the same range. Hypothesis (b) is not ruled out outright (a different, more
specific formulation — e.g. about a symmetry of the operator rather than the visible motion of one
eigenvalue — could still be right), but this specific, cheap, testable prediction it made does not
hold up. (a) and (d)/(e) are relatively more favored by elimination; (a) (integrability-breaking) is
still the best-motivated candidate pending §3.1/§3.5.

### 3.7 Check the other two models (XXZ, tricritical)
Both are already in the "four models" barrier campaign (CLAUDE.md §18) and both are non-integrable
or otherwise structurally distinct from Ising. Does XXZ's transfer matrix (away from the already-
understood exact Néel `P`-sector degeneracy — e.g. using a different, `P`-symmetric quench state if
one exists) or the tricritical model's transfer matrix show the same ≈π-pairing pattern? If the
pairing appears generically whenever integrability is broken (regardless of *which* integrability-
breaking term is responsible), that is further evidence for hypothesis (a); if it is Alcaraz-NNN-
specific, that points back toward something particular to the NNN coupling structure (hypothesis (d)
or a variant of it).

### 3.8 Pending: finish the RDM-vs-RTM check at the wall (T=10) — currently an unbacked claim

A smaller, separate, more immediate item, surfaced while double-checking NB5 rather than while
deriving the pairing candidates above — but worth recording here since it bears on the same
near-degenerate-eigenvector conditioning question.

NB5's "Is the dip a truncation artifact? RTM vs RDM" section (cell 8) asserts that the RDM
truncation, despite converging faster and in fewer iterations than RTM in the pre-wall window
(`T=6,8`), *still* fails once the wall is reached: "a cold start at `T=10` can return a spurious
`|θ0|~2` where the warm-started ladder gives `~1.55`." **Checking `results/data/nb5_rtm_vs_rdm.jld2`
directly shows this is not backed by any cached result** — the file contains only
`(6.0,:rdm),(6.0,:rtm),(8.0,:rdm),(8.0,:rtm)`. The comparison cell's own `temperatures` list already
includes `[6.0, 8.0, 10.0, 11.0, 12.0]`, but `should_regenerate=false`, so `T=10,11,12` were never
computed for either mode — the run was apparently never finished (possibly abandoned for taking too
long), not deliberately excluded. The "~1.55" half of the claim is consistent with `nb8_master.jld2`'s
*warm-started* `T=10` entry (`|λ0|=1.544`), but no cold-started Alcaraz `T=10` run giving `|θ0|~2`
exists anywhere in the repo's cached data — this specific sentence in NB5 currently has no receipts.

**The pending run**: reuse NB5 cell 7's existing code as-is (no new implementation needed) — set
`should_regenerate=true`, restrict `temperatures = [10.0]` first to control cost, same
`Random.seed!(20260627)`, `k=4`, `maxdim=64`, `itermax=1200`, `stuck_after=200`, both `trunc_mode`s.

**What it would tell us**: if RDM *also* comes back `reason=stuck` at `T=10` (matching RTM), that
directly confirms the wall is genuine and truncation-scheme-independent even at the specific point
the notebook's text claims to have checked — strengthening (with an actual receipt this time) the
existing "the wall is physics, not a truncation artifact" conclusion. If RDM instead *converges*
cleanly at `T=10` with a sane `|θ0|`, that would be a materially new finding not currently reflected
anywhere in the project — worth updating NB5's text either way once it actually runs.

**Cost/risk note**: expensive and not to be treated lightly. Iteration counts already grew steeply
between the two cached points (RTM: 321→513 iters, RDM: 190→375 iters, going `T=6→8`); at `T=10`,
right at the documented gap-closing wall for `p=0.1`, either mode could plausibly burn the full
`itermax=1200` without converging, similar to the tricritical `T=4.5` point that ran 7.8h before
being killed (NB11). Should be run **alone** on this 14GB single-Julia-kernel machine, not
concurrently with any other heavy computation (e.g. wait for the NB 5.5 sweep in this same
investigation to finish first).

### 3.9 The Eq.(4) boundary-exponent anomaly ($x_1\approx1.5$ instead of $1/2$) — likely a companion question to the pairing

Surfaced while fixing NB7's Eq.(3) central-charge extraction (§3, above, was about $c$ from
$\lambda_0$'s own phase curvature; this is about a *different* quantity: $x_1$ from the
$\lambda_1-\lambda_0$ phase gap, Eq. 4 of Carignano–Tagliacozzo). At $p=0$, $x_1(p{=}0)=0.498$ exactly
matches the free-BC Ising target and $a_1\approx\pi/4$ as expected. At $p=0.1$, $x_1\approx1.5$ — a
factor of $\sim3$ off, not a small deviation.

Unlike the Eq.(3) $c$-extraction (§ above), this does **not** look like the same kind of selector
artifact: NB7's $\lambda_1$-selector (`lam1_cft`, phase-closest to the physical $\lambda_0$, which by
construction tends to reject the $\approx\pi$-away $-\lambda_0$ partner) was checked directly against
the cache and found to track a **stable** branch across the entire fit window $T=2..11$ — it only
destabilizes at $T=12$, already deep past the documented wall. So whatever is behind the $x_1$ anomaly
survives well inside the region where the eigenvalue bookkeeping is trustworthy, which points toward
it being a **genuine feature of the frustrated model** rather than a bug — but there is currently no
explanation for it.

**Why this belongs next to the pairing question (§1-§2 above).** Eq. (4)'s $\lambda_1$ is, by
definition, the *next* member of the transfer spectrum beyond the leading pair $(\lambda_0,\text{its
partner})$ — exactly the object whose approach and eventual clustering (§1's empirical finding: more
of the $k=4$ leading eigenvalues satisfy the pairing criterion as $T$ grows) is already under
investigation. It is plausible that $x_1$'s anomaly and the pairing mechanism share a common origin
(e.g. hypothesis (a) or (b) in §2), or that $x_1$'s deviation is itself a *diagnostic* for how the
pairing distorts the "next excited state" away from the naive CFT expectation. This is speculative,
not established.

**Concrete next step**: once a candidate mechanism from §2 is identified (most likely via §3.1's
fine-$p$ sweep or §3.5's explicit chiral-operator check), re-examine whether it predicts a *specific*
$x_1$ shift, and check that prediction against this $x_1\approx1.5$ number. Absent that, a cheaper
first move is simply to check whether $x_1(p)$ varies smoothly with $p$ near $p=0$ (reusing §3.1's
planned sweep) — if $x_1$ departs from $1/2$ at the same rate the pairing turns on, that is direct
evidence the two are linked.

### 3.10 The `p=0.3`/`p=0.5` physical-gap anomaly (surfaced while merging NB 5.5 into NB3)

**Observed** (now documented in notebook 3): the physical gap $|\lambda_1|/|\lambda_0|$ for the
Alcaraz model rises smoothly and monotonically with $T$ at `p=0` and `p=0.1`, but is genuinely
erratic at `p=0.3` and `p=0.5` — e.g. at `p=0.3` it goes $0.81\to0.99\to0.86\to0.98\to0.99\to0.99$
over $T=2\ldots7$, jumping up and back down rather than climbing steadily. This uses the same
cheap, exploratory block-PM budget (`itermax=2000`, vs notebook 7's production `itermax=8000`) as
the rest of that cross-check sweep.

**What's confirmed**, from a direct investigation of the cached `k=4` eigenvalues
(`results/data/nb55_pgap.jld2`):

- The `pick_phys`-style selector (restricted to choosing only between the two largest-modulus Ritz
  values, since `block_transfer_eigs` always returns `theta` sorted by modulus) is not innocuous —
  an unrestricted version that searches all 4 candidates for whichever is closest to the previous
  step's value genuinely picks a *different* index at several `T`. But the unrestricted version's own
  "gap" exceeds 1 at some of those points, which is impossible for a genuine dominant/subdominant
  pair — it is tracking a *subdominant* branch, not a validated better choice. Selector fragility is
  real here, but there is no working fix yet.
- The 4-eigenvalue spectral spread is **not** simply "the near-degenerate band forms earlier" for
  larger `p` — the spread `(max−min)/max` is non-monotonic in `T` for `p=0.3` (`0.24→0.08→0.19→
  0.03→0.02→0.03`), and never gets particularly tight at all for `p=0.5` (`0.17→0.17→0.10→0.08→
  0.08→0.09`, staying an order of magnitude looser than `p=0.1`/`p=0.3` ever reach).
- Non-convergence (`reason="stuck"` under the cheap budget) explains the anomaly for `p=0.5`:
  restricting to only the `converged` points there (`T=2,3,5`) recovers a clean monotonic sequence
  (`0.834→0.870→0.927`). It does **not** explain `p=0.3`: the sharp `0.991→0.855` drop occurs
  between two points that are *both* `reason="converged"`.

**Next step**: disentangle the two candidate contributors (selector fragility vs. under-convergence)
by re-running the `p=0.3`/`p=0.5` ladder at notebook 7's production budget (`itermax=8000`,
`stuck_after=400`) and/or a finer `ΔT` (e.g. `0.5` instead of `1.0`) to see whether either change —
alone or combined — restores monotonicity. If neither does, that would point toward a genuine,
not-yet-understood physical effect at stronger frustration rather than a numerical artifact.

---

## 4. Why this matters for the thesis

The barrier campaign (CLAUDE.md §17/§18) already establishes *phenomenologically* that frustration
closes the transfer-matrix gap faster and that the resulting near-degeneracy is what bounds the
method's reach (the "wall"). What is currently missing is a *mechanistic* explanation for why the
near-degeneracy takes the specific form of a tightening ≈π-separated pair, and why this is present
for the non-integrable model but completely absent for the integrable one under the identical
pipeline. If one of the candidates above pans out — most plausibly (a)/(b) given the control
experiment — it would upgrade "the barrier closes faster for frustrated models" from an observed
correlation to a derived consequence of integrability-breaking, which is a stronger and more
citable claim for the thesis's discussion of *why* the barrier appears where it does, not just
*that* it does.
