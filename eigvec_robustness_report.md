# Robust "physical eigenvector" extraction near the transfer-matrix degeneracy
### Research report — candidate approaches, literature, and minimal experiments
*(for discussion. Companion working notebook: `NBs/13_eigvec_robustness.ipynb` (Alcaraz arm, in
progress). Companion code: `block_transfer_eigs`, `src/transverse_tools.jl:142–354`; failed repairs
documented in NB5 §4–5 and CLAUDE.md §17.)*

---

## HANDOFF STATUS (2026-07-12, ~03:05) — read this first if picking up cold

**What exists:**
- This report = the research write-up (candidate approaches A–F, literature, experiment specs E1–E7).
- `NBs/13_eigvec_robustness.ipynb` = the implementation. §0 problem statement, §1 free cache
  diagnostics (no compute), §2 helper functions + the instrumented ladder driver, §3–§6 analysis
  cells for E1/E3/E4/E6/E7 (one cell each), §7 verdicts table (currently empty placeholders — fill
  once the ladder completes). The notebook is the single source of truth; a driver script is
  auto-extracted from its `NB13-HELPERS`/`NB13-DRIVER` marker cells for detached execution (the
  extraction script and the extracted `.jl` are session-scratchpad files, not in the repo — if you
  need to regenerate them, the marker-cell convention is documented in the notebook itself, first
  line of each of those two cells).
- Original plan file (this session only, not repo-tracked):
  `~/.claude/plans/proud-giggling-barto.md` — scope/sequencing rationale, superseded by the
  notebook + this file as the source of truth for content.

**Compute status — a TWO-PHASE run is CURRENTLY QUEUED, not running:**
- Phase 1 (cheap, runs first): the §2b **p-control** — T=3, k=6, p ∈ {0.0, 0.1, 0.3, 0.5}, cache
  `results/data/nb13_pcontrol.jld2` (see the cross-session synthesis in the UPDATE section below
  for why this was added).
- Phase 2: the main ladder — Alcaraz VD2, p=0.1, λ=1, dt=0.1, nbeta=4, k=4, T=2.0:1.0:12.0, same
  budget as NB7's master sweep (itermax=8000, stuck_after=400), warm-started, measuring baseline
  + E1 + E3 + E4 + E6 + E7 in-kernel at every T.
- Cache: `results/data/nb13_eigvec_ladder.jld2` — **currently has only T=2.0** (regression-checked
  against `nb8_master.jld2`, machine-precision match).
- A detached, self-gating launcher (PID 1185282, its own session leader, survives any Claude Code
  session ending) is polling every 2 minutes, waiting for (a) no other `.jl` driver running out of
  any `/tmp/claude-*/scratchpad/` directory, and (b) ≥6.5 GB available RAM, before it starts the
  T=3..12 remainder. As of this writing it's still waiting — `nb3_p01_stress.jl` (a **user-launched**
  production job, unrelated to this investigation) is occupying the machine.
- **Process-safety note for any session, human or Claude, touching this machine**: do not kill
  `nb3_p01_stress.jl` or any other julia process you didn't personally start — this machine has
  14 GB RAM and no headroom for concurrent heavy jobs (an earlier attempt this session caused a
  near-collision, resolved by making the ladder launcher poll for a genuinely idle machine rather
  than a single expected PID). If you want to check progress: `tail
  /tmp/claude-*/*/scratchpad/nb13_ladder.log` won't be visible across sessions (scratchpad is
  session-scoped) — instead just re-read `results/data/nb13_eigvec_ladder.jld2` for cached T's, or
  `ps aux | grep julia` to see whether the gate has opened (a julia process running
  `nb13_driver.jl` means it's actively computing).
- Once T=3..12 land, run NB13 §3–§7 to fill in the verdicts and update this report's UPDATE section
  (below) with the outcomes, then close the loop: report §7's E1/E2/E3/E4/E6/E7 experiment specs
  are now superseded by the actual (parity-aware) implementation in the notebook — treat the
  notebook as authoritative over the original experiment prose where they differ (mainly: E4's
  cluster selection uses boundary-overlap ranking, not modulus ranking — see the UPDATE below).

**Two empirical findings so far** (from cache diagnostics + one validated T=3 smoke point — see
UPDATE below for detail): (1) the wall looks like a near-exceptional-point (bilinear condition
number κ(λ₀) grows 12→1413 over T=2..6), not just two eigenvalues crossing; (2) candidate A's
boundary vector has *exactly zero* overlap with the odd-parity cluster members, which include the
wall-driving −λ₀ partner — so it's automatically immune to that specific contamination, a sharper
mechanism than the report's generic conditioning argument. Both are preliminary (T≤3 data only).

---

## UPDATE (2026-07-12) — preliminary NB13 findings, before the full ladder

Two things surfaced already, from `nb8_master`/`nb7_condnum` (§1) and a validated T=3 smoke run:

1. **The wall is a near-EP (self-orthogonality), not a plain eigenvalue collision.** The bilinear
   eigenvalue condition number κ(λ₀)=‖L₀‖‖R₀‖/|⟨L₀|R₀⟩| grows 12→43→124→417→**1413** over T=2..6
   (phase rigidity r₀=1/κ → 0.02 already at T=3, uniform across all four members). The leading pair
   is becoming self-orthogonal in the bilinear form — the E1 diagnostic scenario, and plausibly the
   *emergent dual-unitarity* fingerprint itself. This refines §1.2: we are heading toward an
   exceptional point, and the "gap → 1e-3" language should be read as κ → ∞, not just |λ₁|→|λ₀|.
2. **Candidate A gets a much sharper mechanism than the generic divided-difference argument.** The
   structured `fw_tMPS` boundary (|X⁺⟩, even parity) has boundary overlaps ⟨Lⱼ|b⟩ = **[3.11, 0, 0,
   1.21]** — *exactly zero* on the odd-sector members {2,3}, which include the −λ₀ partner whose
   modulus-degeneracy with λ₀ drives the wall. So E^ℓ·b is automatically **parity-projected onto the
   even sector**, structurally excluding the contaminating partner. The prediction sharpens: the wall
   is even↔odd mixing of the *individual* Ritz vector during iteration; a boundary-projected object
   is immune by symmetry, not just by conditioning. NB13's candidate-A cluster is now selected by
   boundary overlap (`:bpair`) rather than by modulus. Whether the two *even* boundary-visible members
   {1,4} themselves become modulus-degenerate at the wall (making ℓ genuinely matter) is the open
   question the full ladder resolves.

**3. Cross-session synthesis (later on 2026-07-12) — the p=0.1 anomaly and the tower/partner
classification.** A parallel session (NB3 gap-anomaly investigation) found p=0.1 is *spectrally
anomalous*: by exact dense diagonalization, the x≈0.5 boundary descendant (C-T Eq.4's free-BC
state) is amplitude-suppressed ~5× at p=0.1 — a sharp **non-monotone** dip (normal at p=0, 0.3,
0.5), not a k-truncation artifact (confirmed at k=8/10, `nb3_p01_stress.jld2`). They classify the
block spectrum by phase: ph(z)=angle(−z), Δφ to λ₀ wrapped to (−π,π], x=(vT/π)|Δφ|; |Δφ|<π/2 =
"tower" (λ₀ + CFT descendants), else "partner" (π-shifted near-equal-modulus copies). Reclassifying
**my** T=3 smoke block in their convention interlocks the two pictures exactly:

  | member | \|θ\|/\|θ₀\| | Δφ/π | x | \|⟨L\|b⟩\| | class |
  |---|---|---|---|---|---|
  | 1 (λ₀) | 1.000 | 0 | 0 | 3.11 | tower |
  | 2 | 0.978 | −0.94 | — | **0 (exact)** | partner |
  | 3 | 0.892 | −0.81 | — | **0 (exact)** | partner |
  | 4 | 0.864 | +0.25 | **1.501** | 1.21 | tower |

  So: (a) their "partner" class *is* my hard-zero-boundary-overlap sector — two independent
  definitions, same split; (b) the boundary couples to λ₀ + the **x≈1.5** state because the x≈0.5
  state is absent from the k=4 block at p=0.1 (their suppression) — which is very likely the
  long-standing NB7 x₁≈1.5 anomaly, now visible in boundary-overlap language; (c) their hard/soft
  distinction is worth preserving: the **hard zero** on partners is structural (parity — why
  candidate A is immune to partner contamination at any p), while the **soft ~5× dip** of the x≈0.5
  tower amplitude is p≈0.1-specific and lives *within* the boundary-visible sector (it does not
  threaten candidate A's mechanism, but it makes p=0.1 a non-generic anchor for boundary-overlap
  geometry). Consequences drawn: NB13 gained a **§2b p-control** (T=3, k=6, p ∈ {0.0, 0.1, 0.3,
  0.5}; cache `nb13_pcontrol.jld2`; queued *before* the main ladder in the detached driver) that
  measures boundary overlaps, per-member rigidity, the phase/x classification, and the E6 gradings
  across p — testing three falsifiable predictions: partners keep exact-zero boundary overlap at
  every p; the partner class vanishes at p=0 (their finding — also *their* prediction for my E6:
  the grading coupling should be ~0 at p=0 and grow with p, a free test of the D-staggering
  hypothesis); the x≈0.5 tower amplitude dips only at p=0.1. The p=0.1 ladder stays the headline
  (it is where the wall and the thesis result live), but its boundary-overlap interpretation is now
  explicitly conditioned on this control.

Everything below is the original research report; the points above are the first empirical
contact with it.

---

## 0. Executive summary

**The one structural fact everything follows from:** near a cluster of eigenvalues with internal
gap `g`, an *individual eigenvector* has condition number ~ `ε/g` under a perturbation `ε` of the
operator — but (i) the *invariant subspace* of the whole cluster is conditioned by the separation
`sep` of the cluster from the *rest* of the spectrum (O(1) here), and (ii) **any analytic matrix
function of E applied to a physical vector — e.g. `E^L·b` — is conditioned by the derivative of the
function, with no `1/g` factor at all** (divided-difference bound, Higham, *Functions of Matrices*,
2008). The eigenvalue route survives the wall because eigenvalues of a diagonalizable-but-clustered
matrix are well-conditioned; the entropy route dies because it asks the one question — "which
direction inside the cluster?" — that is mathematically ill-posed at the degeneracy. Any proposal
must be graded against this rule: does it stop asking that question, or does it just re-ask it in
different clothes?

**Ranked verdicts** (details in the numbered sections):

| # | Approach | Verdict | Cost | Section |
|---|----------|---------|------|---------|
| A | **Finite-L boundary-weighted cluster combination** `R_phys(L) ∝ Σᵢ θᵢᴸ ⟨Lᵢ\|b⟩ Rᵢ` — a *pure* state, gauge-invariant over the cluster, = the physical RTM of a width-L system | **Sound. Top candidate.** Replaces the ill-posed T→∞ single-vector question with a well-posed matrix-function evaluation; reduces exactly to (L₀,R₀) when `g·L ≫ 1`; smooth through the wall; turns the wall into a physical crossover-scaling result | Post-processing only: k overlaps + one `lincomb_mps` + `gen_renyi2` per (T,L) | §5.1 |
| B | **Symmetry-resolved cluster rotation** — diagonalize a commuting symmetry Q (not E) inside the cluster; Q's cluster eigenvalues are O(1)-separated → well-conditioned | **Sound where a symmetry drives the degeneracy** (XXZ Z₂ definitely; Alcaraz ±pair if the staggering hypothesis holds). Needs one derivation: the symmetry's action on the *temporal* legs | k² cheap overlaps + a k×k diagonalization, once per T | §5.2 |
| C | **Refined Ritz extraction** with the known θ_phys (Jia): minimize ‖(E−θ)x‖ over the block instead of reading eigenvectors off `eigen(W)` | **Cheap hygiene, cannot beat 1/g.** Fixes exactly the symptom class "values converge, vectors don't"; removes the `pinv(S)`/`eigen(W)`/`pinv(VR)` fragilities; improves constants only. Also a *diagnostic*: if the wall moves, part of it was algebra-ε | ~3k² extra Hermitian `inner()` calls, **final iteration only** | §2.2 |
| D | **Subspace-Procrustes ("diabatic") T-ladder tracking** — align the whole current cluster *basis* to the previous T's basis (k×k rotation), not one stale vector into a contaminated basis | Mechanistically fixes failed attempt #2's failure mode; but the smooth frame is smooth *by construction* — path-dependence (Mead–Truhlar) means no guarantee it carries the *right* vector out the other side. Use as a frame for A/B, not standalone | k² overlaps per T-step | §5.3 |
| E | **Inverse iteration / RQI seeded with θ_phys** (full MPS space) | **Moves the ill-conditioning; not worth it.** Converges beautifully — to the exact eigenvector of the *truncation-perturbed* operator, which differs from truth by the same ε/g. Also expensive (MPO linear solves = DMRG-sweep-scale) | High | §2.1 |
| F | **Hermitian dilation [[0,M],[M†,0]]** (question 2b) | **Sound math, wrong physics.** Dilating E gives singular vectors = the *folded/Hermitian* temporal picture, whose entropy grows linearly in T — the exact object the RTM method exists to avoid. Dilating (E−θ) recovers the eigenvector but re-imports the same 1/g through the singular-value gap. Useful as a *proof* that the ill-conditioning is intrinsic to the question | — | §3 |

**Recommended sequence:** run the diagnostics (E1, E2 in §7) to establish how much of the wall is
implementation-ε vs physics; implement A (it is pure post-processing on already-converged blocks);
pursue B for XXZ where the symmetry is explicit; keep C as a cheap upgrade to the final extraction.

---

## 1. Background: the conditioning structure of the problem

### 1.1 Individual eigenvectors vs invariant subspaces vs matrix functions

For a diagonalizable non-Hermitian E with simple eigenvalue λ₀, first-order perturbation of the
right eigenvector under E → E+εF is

  δR₀ = Σ_{j≠0} (Lⱼᵀ F R₀)/(λ₀−λⱼ) · Rⱼ  →  the intra-cluster term carries the `1/g`.

Three objects with progressively better conditioning (Stewart & Sun, *Matrix Perturbation Theory*,
1990; Stewart, SIAM Rev. 1973; Kato 1966):

1. **Individual eigenvector**: κ ~ 1/g. This is what the current Ritz-pair entropy uses. Doomed at
   the wall by mathematics, as already established by the exact-diagonalization ground truth (NB5).
2. **Cluster invariant subspace / spectral projector**: κ ~ 1/sep(cluster, rest). O(1) here — the
   4-fold cluster stays well-separated from θ₅ downward. This is why the *span* of the computed
   block stays good even when its individual members scramble.
3. **Analytic matrix functions of E**: for f analytic on a neighborhood of the spectrum,
   ‖f(E+εF) − f(E)‖ is controlled by divided differences of f, e.g. for f(z)=z^L the intra-cluster
   divided difference is (θ₀ᴸ−θ₁ᴸ)/(θ₀−θ₁) ~ L·θᴸ⁻¹ — **bounded, no 1/g** (Higham 2008, Ch. 3).
   The eigendecomposition of f(E) is as ill-conditioned as E's; but *evaluating f(E) on a vector*
   never passes through the eigenbasis.

The design rule: **an observable is safe iff it is expressible as a gauge-invariant function of
(E, cluster subspace, external physical data) without resolving individual intra-cluster
directions.** Graded against this rule:

- *Failed attempt #1* (projector P = Σ|Rᵢ⟩⟨Lᵢ|): gauge-invariant ✓, well-conditioned ✓ — but it is
  f(E) = 1 on the cluster, i.e. the *infinite-temperature* state of the cluster. Its entropy adds
  log(m) classical mixing by construction. Right conditioning class, wrong function.
- *Failed attempt #2* (project previous vector onto current block): asks the individual-direction
  question with inputs that are already scrambled, and anchors on a single stale vector. Both the
  question and the mechanism were wrong.
- *Candidate A below*: f(E) = E^L applied to a physical boundary vector — right conditioning class
  *and* the physically correct weighting. It is attempt #1 with the correct f.

### 1.2 Diagnostic: avoided crossing vs exceptional point (do this first)

Non-normality adds a second failure axis: at an **exceptional point** (EP) the eigenvectors
*coalesce* (⟨L₀|R₀⟩_noconj → 0), the matrix becomes defective, perturbation theory switches from
Taylor to Puiseux (√ε) scaling (Kato 1966; Moro–Burke–Overton, SIMAX 18, 1997), and even the 2×2
cluster restriction is non-diagonalizable — candidates A/B still work (they never diagonalize the
cluster if implemented via the restriction matrix, see §5.1 caveat), but any eigenvector-based
language breaks down entirely.

Your data argues *against* a true EP: eigenvalues correct to 1e-8..1e-13 through the wall means the
eigenvalue condition number κ(λ) = ‖L‖‖R‖/|⟨L|R⟩_noconj| stays modest (it diverges at an EP). But
this should be *measured*, not inferred: the **phase rigidity** rⱼ = |⟨Lⱼ|Rⱼ⟩_noconj|/(‖Lⱼ‖‖Rⱼ‖)
(Rotter, J. Phys. A 42, 153001, 2009) per cluster member across the T-ladder is a one-line
diagnostic on data you already have, and it discriminates the two scenarios:

- r stays O(1) through the wall → clustered-but-diagonalizable; A/B/C all apply as stated.
- r → 0 approaching the wall → near-EP; the ± pair is trying to merge into a Jordan block; the
  final bi-orthonormalization `(1/√⟨L|R⟩)·` (transverse_tools.jl:344) is itself blowing up, and the
  cluster restriction must be handled in Schur form, never eigen form.

NB7 already computes κ(λ₀)=‖L₀‖‖R₀‖ for the sym-vs-asym comparison — this is the same quantity;
it just needs to be read per-member across the wall.

### 1.3 Where ε lives in the current implementation

The eigenvalues are good to 1e-8 while the effective vector-level contamination is ~1e-3 (O(1)
rotation at g~1e-3). That mismatch localizes the dominant ε in the MPS-level operations, not the
k×k algebra:

- `applyn` truncation of each MPO application (cutoff-scheduled);
- the joint RTM `truncate_sweep` — your own comment at transverse_tools.jl:249: "the non-Hermitian
  SVD is ill-conditioned right at the gap closing" — this injects noise *precisely into the
  rotating intra-cluster directions*, every iteration;
- `lincomb_mps` de-mixing each iteration (mitigated by `basis=:schur`, whose rotations are unitary);
- the k×k algebra (`pinv(S)` at rtol 1e-12, `eigen(W)`, `pinv(VR)`) — likely subdominant, but
  candidate C removes it for free.

This matters because the wall's *position* is set by where ε/g ~ 1. If ε has a reducible
implementation component, the wall moves without any new theory. Experiment E2 (§7) quantifies this.

---

## 2. Question 1: extracting the vector given θ_phys

### 2.1 Inverse iteration / RQI in the full MPS space — verdict: moves the ill-conditioning

The idea is seductive: θ_phys is known to ~1e-8 while g ~ 1e-3, so one inverse-iteration step with
shift θ_phys amplifies the target eigenvector over its cluster partner by |θ−λ₁|/|θ−λ₀| ~
g/δθ ~ 10⁵. Two fatal problems:

1. **It converges to the eigenvector of the wrong operator.** The operator actually available is
   Ẽ = E + (truncation), with ‖Ẽ−E‖ ~ ε_trunc. Inverse iteration finds Ẽ's eigenvector *exactly*
   — which differs from E's by the same ε_trunc/g rotation. Knowing θ to 1e-8 reduces the *solver's*
   contribution to the error; it cannot reduce the vector's sensitivity to perturbations already
   baked into the operator. (This is the general pattern for this whole question: the eigenvalue's
   precision is not a resource that transfers to the eigenvector, because their condition numbers
   are governed by different quantities — |⟨L|R⟩| vs 1/g.)
2. **Cost**: (E−θ)⁻¹ applied to a tMPS is an MPO linear solve — a DMRG-style sweep with a shifted
   operator, i.e. at least "a second full block sweep" per iteration, your stated non-starter.

For the bilinear structure specifically: inverse iteration is well-defined (shift-invert of a
general matrix), and the left vector needs the transposed solve (E−θ)⁻ᵀ — consistent with your
`mpoT = swapprime(mpo,0,1)` convention. No conceptual obstruction, just no payoff.

### 2.2 Refined Ritz extraction inside the existing block — verdict: cheap, do it, expect constants

This is the legitimate version of "use the known θ": **refined projection** (Jia, *Linear Algebra
Appl.* 259, 1997; Jia, *Comput. Math. Appl.*, 2001 — see the [residuals of refined projection
methods paper](https://www.sciencedirect.com/science/article/pii/S0898122100003217) and the
[convergence analysis of Ritz vs refined Ritz vectors](https://www.researchgate.net/publication/2450993_On_the_Convergence_of_Ritz_Values_Ritz_Vectors_and_Refined_Ritz_Vectors)).
The classical result matches your symptom exactly: **for non-Hermitian problems, Ritz values can
converge while Ritz vectors fail to converge or are non-unique; the refined Ritz vector — the
minimizer of ‖(E−θI)x‖ over the subspace at the (accurately known) θ — converges unconditionally
as the subspace converges.** For a *cluster* of nearby Ritz values, the refined vectors all converge
to the cluster eigenspace (recent treatment in the
[refined CJ–SS–RR method paper, arXiv:2605.12846](https://arxiv.org/pdf/2605.12846)), which is
the honest statement of the limit: refinement gets you cleanly *into* the right subspace; it cannot
pick the direction *within* it any better than the gap allows.

Concrete recipe in your setting (right vector; left is mirror with ATL and Eᵀ):

- Basis: the current applied block AR₁..AR_k (already computed). Seek x = Σⱼ cⱼ Rⱼ minimizing
  ‖E x − θ_phys x‖₂ / ‖x‖₂.
- This needs three k×k **Hermitian** Gram matrices (genuine `inner()`, with conjugation — this is a
  2-norm minimization, deliberately outside the bilinear structure):
  G = [⟨Rᵢ|Rⱼ⟩], X = [⟨Rᵢ|ARⱼ⟩], H = [⟨ARᵢ|ARⱼ⟩].
- Minimize c†(H − θ̄X − θX† + |θ|²G)c subject to c†Gc = 1: smallest eigenpair of a k×k Hermitian
  PSD pencil. ~3k² = 48 extra MPS inner products at k=4, **needed only once, after convergence**,
  not per iteration.
- Bonus: it bypasses `pinv(S)`, `eigen(W)`, and `pinv(VR)` entirely for the final extraction, and
  gives a residual norm ‖(E−θ)x‖ as a free, honest error bar on the vector.

Expected outcome, stated for the record so the experiment is falsifiable: **if the wall is purely
physics (ε dominated by truncation), the refined-Ritz entropy breaks at the same T; if part of the
wall was k×k algebra noise, it moves outward by a bit.** Either result is worth having.

---

## 3. Question 2b: the Hermitian dilation — sound math, wrong physics

The proposal: embed into H = [[0, M],[M†, 0]], whose eigenpairs are (±σᵢ, [uᵢ; ±vᵢ]) — singular
triplets of M — and which is Hermitian, hence has well-conditioned eigenvectors whenever the
singular values are separated. Two ways to deploy it, both assessed:

**(i) Dilate E itself.** Correct as linear algebra, but the singular vectors of a *non-normal* E
are unrelated to its eigenvectors (they coincide only for normal matrices — and the transfer
matrix's non-normality is the whole game here). Physically it is worse than unrelated: the
SVD/Hermitian structure of the temporal wavefunction is exactly the **folded picture**, and the
entropy it defines is the standard Hermitian temporal entanglement entropy, which **grows linearly
in T after a quench** — this is precisely the complexity result of Carignano–Marimón–Tagliacozzo
(PRR 6, 033021, 2024) that motivated the RTM/generalized-entropy construction in the first place.
The dilation's eigenvectors are well-conditioned *because they describe different physics*, and
that physics is the computationally hostile one. Additionally, the dilation replaces the bilinear
pairing (`overlap_noconj`, transpose) with the sesquilinear one (adjoint) — the transition-matrix
object |R₀⟩⟨L₀| is bilinear by construction, so nothing built from M† can represent it directly.

**(ii) Dilate the *shifted* operator (E − θ_phys·I).** This one is genuinely clever and worth
recording why it still fails: the null right singular vector of (E−θ₀) *is* the right eigenvector,
and the null left singular vector is the conjugate of the transpose-left eigenvector (E ᵀu = θu ⇔
(E−θ)ᵀu = 0), so one SVD-type problem would deliver the *pair*, with Hermitian conditioning. But
Hermitian conditioning of a singular vector is governed by the gap between *its* singular value and
the next one — and σ₂(E−θ₀I) − σ₁(E−θ₀I) ≈ |λ₁−λ₀| = g (up to non-normality factors). **The 1/g
re-enters through the singular-value gap.** This is the cleanest available proof that the
ill-conditioning is intrinsic to the question "which direction inside the cluster", not to any
particular algebraic formulation of it. Recommend keeping this argument in the thesis's methods
discussion; do not implement.

---

## 4. Question 2: literature survey — how other fields track through near-degeneracies

### 4.1 Numerical linear algebra

- **Invariant-subspace perturbation, `sep`, clusters**: Stewart (SIAM Rev. 15, 1973), Stewart & Sun
  (*Matrix Perturbation Theory*, 1990). The canonical statement of §1.1. The prescribed object for
  a cluster is the **ordered Schur / deflating subspace**, never the eigenvector basis — your
  `basis=:schur` option is exactly this move, applied to the iteration; the proposals here extend
  it to the *extraction*.
- **Defective/nearly-defective perturbation theory**: Moro, Burke, Overton (SIAM J. Matrix Anal.
  18, 1997) — Lidskii–Puiseux perturbation of Jordan structure; Wilkinson's distance to the nearest
  defective matrix. Relevant if the phase-rigidity diagnostic (§1.2) says "near-EP".
- **Smooth continuation of eigendecompositions along parameter paths**: Dieci & Eirola, "On smooth
  decompositions of matrices" (SIAM J. Matrix Anal. 20, 1999); Dieci & Friedman, "Continuation of
  invariant subspaces" (Numer. Lin. Alg. Appl. 8, 2001); see also
  [continuous decompositions and coalescing eigenvalues for parameter-dependent matrices](https://www.researchgate.net/publication/282766494_Continuous_Decompositions_and_Coalescing_Eigenvalues_for_Matrices_Depending_on_Parameters).
  Core message: along a one-parameter path, what can be continued smoothly through a
  near-degeneracy is the **subspace**, with a smooth gauge chosen by an ODE/Procrustes condition —
  individual eigenvectors cannot. This is `pick_phys`'s philosophy lifted from eigenvalues to
  subspaces, and is the theoretical backing for candidate D.
- **Very close recent analog**: predictor–corrector homotopy continuation of *non-Hermitian*
  dispersion curves through EP-infested regions, using an auxiliary Hermitian "anchor" problem to
  assign global branch identities ([arXiv:2605.15089](https://arxiv.org/pdf/2605.15089)). Their
  finding that "local trackers fail silently when EPs approach the tracking path" is a useful
  caution for `pick_phys` too: worth checking that no EP sits close to the T-ladder path (the
  phase-rigidity scan does this).
- **Refined and harmonic Ritz vectors**: §2.2 references; also
  [harmonic Ritz convergence analysis](https://arxiv.org/pdf/1603.01785). Harmonic Ritz targets
  interior eigenvalues — not needed here (the cluster is exterior/dominant); refined Ritz is the
  relevant one.

### 4.2 Quantum chemistry: the adiabatic→diabatic playbook

At an avoided crossing, adiabatic eigenvectors rotate violently (the same ε/g); a century of
practice says: **stop diagonalizing, rotate the 2-state cluster to a basis defined by a smooth
auxiliary criterion.** Key facts with citations:

- Strictly diabatic bases **do not exist** in general (Mead & Truhlar, J. Chem. Phys. 77, 6090,
  1982) — the derivative coupling has a nonremovable part; only *quasi*-diabatization is possible.
  Translation to our problem: no exact smooth-frame construction exists either; any smooth frame
  carries path-dependent error. This is the fundamental caveat on candidate D.
- Quasi-diabatization strategies, in decreasing relevance to us:
  (a) **Property-based**: diagonalize a physical *property operator* (dipole, charge) within the
  near-degenerate block, because its eigenvalues stay O(1)-separated where the energy's don't
  (Werner & Meyer 1981; block-diagonalization of Pacher–Cederbaum–Köppel, J. Chem. Phys. 89, 1988;
  review: Van Voorhis et al., Annu. Rev. Phys. Chem. 61, 149, 2010). **This is exactly candidate B**
  with the property operator = the Z₂ symmetry / staggering operator.
  (b) **Overlap/Procrustes-based**: maximize overlap of the current block with the previous
  geometry's block (Löwdin orthogonalization of the mixed overlap matrix). **This is candidate D**,
  and the chemistry experience is directly transferable: works well through a single crossing,
  accumulates gauge drift over long paths.

### 4.3 Non-Hermitian / exceptional-point physics

- Kato (1966) — the mathematical foundation; Heiss, J. Phys. A 45, 444016 (2012) — the physics of
  EPs; Ashida, Gong, Ueda, Adv. Phys. 69 (2020) — modern review; Rotter (J. Phys. A 42, 153001,
  2009) — **phase rigidity** as the standard order parameter for approach to an EP (§1.2
  diagnostic). General background:
  [Exceptional points in non-Hermitian systems](https://www.emergentmind.com/topics/exceptional-points-eps),
  [perturbation theory in the complex plane](https://iopscience.iop.org/article/10.1088/1361-648X/abe795).
- The EP literature's blunt lesson: at/near an EP the individual-eigenvector question has *no*
  answer (the Riemann-sheet structure means the "same" eigenvector returns as the *other* one after
  encircling); everything well-defined is a function of the 2-dim Jordan pair. Consistent with §1.1.

### 4.4 Tensor networks: degenerate transfer-matrix fixed points

- **Transfer matrices and excitations with MPS**: Zauner-Stauber, Haegeman, Verstraete et al.,
  [NJP 17, 053002 (2015)](https://iopscience.iop.org/article/10.1088/1367-2630/17/5/053002) — the
  transfer-matrix spectrum ↔ correlation structure dictionary.
- **Symmetry-broken phases = exactly degenerate dominant transfer eigenvalues**: standard iDMRG/
  VUMPS practice is *never* to fight the degeneracy numerically — one works with
  **symmetry-resolved (sector-labeled) transfer matrices** or picks the minimally-entangled,
  maximally-symmetry-broken states; with spin-flip-symmetric states the eigenvalues come in exact
  degenerate pairs and one either imposes the symmetry explicitly or selects the sector by a
  reference/boundary state (see e.g.
  [Z_N symmetry breaking in PEPS](https://arxiv.org/pdf/1703.04137),
  [long-range order and symmetry breaking in PEPS](https://arxiv.org/pdf/1505.04217),
  [iDMRG limit cycles / infinite boundary conditions](https://arxiv.org/pdf/1804.09163)).
  This is field-internal precedent for **candidate B** (sector labels) and **candidate A**
  (boundary-state selection of the physical combination) — your XXZ 4-fold cluster (2 sectors × ±)
  is structurally identical to a symmetry-broken 2D-classical boundary-MPS problem.
- The iDMRG "limit cycle" literature is also the likely home of your open −λ₀-partner question:
  eigenvalue pairs (λ, −λ) of a one-column transfer matrix are the fingerprint of **period-2 /
  antiperiodic structure** (band folding — a two-column blocking maps both onto λ²). If the ±
  partners are related by a temporal staggering operator D (diagonal ±1 pattern on the time legs),
  then D is a cluster-resolving symmetry for candidate B. Probe: compute |⟨L₁|D|R₀⟩| — O(1) confirms
  the hypothesis (experiment E6).

---

## 5. Question 3: alternative observables

### 5.1 Candidate A — the finite-L, boundary-weighted pure-state combination (top pick)

**Definition.** After `block_transfer_eigs` converges (cluster members bi-normalized, ⟨Lᵢ|Rᵢ⟩=1),
choose a physical boundary tMPS `b` and a spatial width L, and form the **pure** pair

  |R_phys(L)⟩ ∝ Σᵢ (θᵢ/θ₁)ᴸ ⟨Lᵢ|b⟩_noconj |Rᵢ⟩  (sum over the cluster, m=2 or 4 terms)
  ⟨L_phys(L)| ∝ Σᵢ (θᵢ/θ₁)ᴸ ⟨b'|Rᵢ⟩_noconj ⟨Lᵢ|

then S₂(T; L) = `gen_renyi2(L_phys, R_phys)` as usual (the overlap normalization handles the rest).
The (θᵢ/θ₁)ᴸ rescaling avoids overflow; use the high-precision tracked θ's.

**Why it is well-conditioned.** The combination equals P_cluster·Eᴸ·|b⟩ up to normalization — an
*analytic matrix function of E evaluated on a vector*, restricted to the well-conditioned cluster
subspace. It is invariant under any re-mixing of the computed bi-orthogonal cluster basis (the
gauge freedom that scrambles individual Ritz pairs cancels between the ⟨Lᵢ|b⟩ coefficients and the
|Rᵢ⟩ vectors — they transform contragradiently). The only surviving errors are (i) out-of-cluster
leakage of the computed subspace, ~ ε/sep = O(ε), and (ii) θ-errors, ~1e-8·L. **No 1/g anywhere.**
By §1.1 item 3, this is the unique conditioning class that is both safe and physically weighted.

**Why it is the right physics.** For a *finite* spatial chain of width ~L, the transverse
contraction is literally ⟨boundary|Eᴸ|boundary⟩-structured, and the temporal RTM at a cut is built
from exactly this propagated boundary vector — a pure state, not a mixture. The thermodynamic-limit
single-eigenvector object is the L→∞ limit of this, and that limit *does not commute with the
degeneracy*: as g→0 the required L to converge to a single eigenvector diverges as 1/g (this is
your own §17 observation that "a direct finite-L forward contraction needs L~1e4 at the gap
closing" — restated positively: the finite-L object is fine, it is only the limit that breaks).
So:

- **g·L ≫ 1**: the sum is dominated by the tracked θ_phys term → reduces *exactly* to the current
  single-pair entropy. (Provable and testable — experiment E4a.)
- **g·L ≲ 1**: a genuine coherent combination — which is what a width-L system physically has at
  that T. The "wall" becomes a **crossover** at g(T)·L ~ 1, and S₂(T; L) is a smooth two-parameter
  family. If the CFT dome is universal, the family should collapse in the scaling variable g·L —
  a potentially thesis-grade result ("the entanglement barrier as a finite-size crossover"), rather
  than a workaround.

**How it differs from the two failed attempts** (per your ground rule):
- vs. attempt #1 (projector): the projector is the *equal-weight, incoherent* (mixed, rank-m)
  cluster object — f(E)=1; this is the *physically-weighted, coherent* (pure, rank-1) one —
  f(E)=Eᴸ·|b⟩. No classical mixing entropy: it's one MPS, not a sum of dyads.
- vs. attempt #2 (continuity projection): no stale previous-T vector anywhere; the weights are
  recomputed at each T from current data by a formula; and the ill-conditioned individual vectors
  never appear alone — only in the gauge-cancelling combination.

**The two open design choices** (flagging honestly):
1. **What is `b`?** The physically motivated choice is the actual spatial-edge column of the
   network. Note the structured `fw_tMPS` seed — the one the drivers deliberately overwrite because
   it is "symmetry-special and gets trapped in a subdominant Z₂ sector" — is plausibly *exactly*
   this object: its sector structure is not a bug here but the physical information about which
   cluster combination a real finite chain populates. Worth checking its cluster overlaps ⟨Lᵢ|b⟩
   first (experiment E4b). Fallback: any generic b with all-nonzero cluster overlaps still gives a
   valid finite-L regularization; the L-dependence study reveals how much the choice matters.
2. **L is a new parameter.** That is a feature, not a bug (it restores a limit that was being taken
   implicitly and illegitimately), but the thesis narrative must present S₂(T; L) as an L-family
   with the g·L crossover, not as "the" entropy.

**Near-EP caveat**: if the phase-rigidity scan says the ± pair approaches an EP, the bi-normalized
eigenbasis expansion above degrades (⟨Lᵢ|Rᵢ⟩→0 before normalization). The EP-safe formulation of
the same object: build the m×m cluster restriction W_c and boundary coefficient vector in the
*Schur* (orthonormal-column) cluster basis, compute W_cᴸ·(coeffs) by repeated m×m multiplication,
then `lincomb_mps` with the resulting coefficients. Identical result when diagonalizable, and never
inverts an eigenvector matrix.

**Cost**: m `overlap_noconj` calls + one `lincomb_mps` + one `gen_renyi2` per (T, L) — pure
post-processing; the L-scan reuses the same converged block. Zero change to the iteration.

### 5.2 Candidate B — resolve the cluster with a symmetry, not with E

**Idea** (= property-based diabatization, §4.2; = sector-resolved transfer matrices, §4.4): when
the near-degeneracy is symmetry-driven, there exists an operator Q with [E, Q]=0 whose eigenvalues
*within the cluster* are separated by O(1) (e.g. ±1 for Z₂). Diagonalizing the k×k matrix
Q̃ᵢⱼ = ⟨Lᵢ|Q|Rⱼ⟩ restricted to the cluster is then a **well-conditioned** rotation (condition set
by Q's O(1) cluster gap, not E's closing one), and it labels the cluster members by sector —
recovering per-sector Ritz pairs whose remaining intra-sector gap is the *slow* barrier
(your |θ₅/θ₁| finding), buying back the reach the symmetry cluster destroyed.

**What Q is in the temporal picture** — the one derivation this needs: a global on-site symmetry
u^⊗N of the spatial propagator (XXZ: the Z₂ spin flip, in the rotated Néel frame still an on-site
product) acts on the spatial MPO's *virtual* bond through a gauge matrix V — by the fundamental MPO
symmetry theorem, u-conjugation of the physical legs equals V-conjugation of the virtual legs
(Cirac–Pérez-García–Schuch–Verstraete MPS/MPO representation theory). After the 90° rotation the
spatial virtual legs *are* the temporal physical legs, so **Q acts on the tMPS as a product of the
same single-site matrix V on every time leg** — trivially cheap to apply. V is found numerically
once per model: solve the one-site fixed-point equation W·(u⊗u†-conjugation) = (V⁻¹⊗V)-conjugation
of W for the given U(dt) tensor (a small linear problem), or read it off analytically for the Murg
construction. Then Q̃ᵢⱼ costs k² overlaps with one on-site operator inserted.

**Assessment.** (a) Fundamentally sound *where a commuting symmetry with O(1) cluster separation
exists*: XXZ's 2-sector structure, yes by construction; the ±λ₀ partners, *if* the staggering
hypothesis (§4.4, E6) holds; the tricritical thickening band, probably not (no known symmetry —
the band is charge/operator-content-driven, and B does not apply). (b) Cost: negligible per T.
(c) Note B alone does not finish the job for XXZ: the physical Néel boundary populates *both*
sectors, so after sector resolution one still needs the physical combination — which is candidate A
applied *within* the symmetry-refined basis. B and A compose naturally: B makes the cluster
block-diagonal with well-separated blocks; A weights the blocks physically.

### 5.3 Candidate D — subspace-Procrustes ("diabatic frame") T-ladder tracking

For completeness, the corrected version of failed attempt #2: at each T-step, compute the m×m mixed
overlap Oᵢⱼ = ⟨Lᵢ^(T-ΔT)|Rⱼ^(T)⟩_noconj between the previous and current *cluster bases* (not one
vector), and rotate the current basis by the orthogonal-Procrustes factor of O (polar
decomposition / Löwdin). Subspace-to-subspace alignment is conditioned by sep, not g — so this
fixes the *mechanism* of attempt #2 (single stale vector into a scrambled basis). The physical
vector is then selected once, in the clean window, and carried by the smooth frame.

Honest assessment: the frame is smooth *by construction*, so smoothness of the resulting entropy is
not evidence of correctness; Mead–Truhlar path-dependence means gauge error accumulates along the
ladder with no internal error signal. Use it (i) as a stress test — its drift vs attempt #2's
0.26→0.34→0.86 quantifies how much of that failure was mechanism vs fundamentals — and (ii) as the
smooth gauge inside A/B, not as a standalone answer.

### 5.4 What cannot exist (worth one paragraph in the thesis)

A single-number, gauge-invariant entropy that (a) resolves the cluster to one state and (b) is
continuous through g→0 **cannot exist without extra data**: the spectral projector of an individual
eigenvalue is a discontinuous (pole-carrying) function of the operator at the degeneracy, so any
resolver must import information beyond E itself. The available imports are exactly: a length scale
(candidate A's L), a symmetry label (candidate B's Q), or a path/history (candidate D / `pick_phys`).
This classifies the solution space and explains *a priori* why attempts #1–2 — which imported
nothing (equal weights) or the wrong thing (a stale vector) — had to fail. Framed this way, the
existing eigenvalue-continuity trick, A, B, and D are the four instances of one principle:
**add the physical selector that the mathematics demands.**

---

## 6. Assessment against your criteria (question 4, consolidated)

| Approach | (a) Fundamentally sound? | (b) Cost inside block PM | (c) Experiment |
|---|---|---|---|
| A: finite-L combination | Yes — matrix-function conditioning, no 1/g; physical for finite width | Post-processing only (m overlaps + lincomb + gen_renyi2 per (T,L)) | E4 |
| B: symmetry resolution | Yes, where symmetry exists (XXZ ✓, ± pair ?, tricritical ✗) | One-time V derivation; then k² overlaps per T | E5, E6 |
| C: refined Ritz @ θ_phys | Improves constants + removes algebra fragility; cannot beat ε_trunc/g | ~3k² Hermitian inners, final iteration only | E3 |
| D: subspace Procrustes | Mechanism sound; gauge drift unbounded in principle | m² overlaps per T-step | E7 |
| E: inverse iteration (full space) | No — converges to perturbed operator's vector; same ε/g | MPO linear solves — prohibitive | skip |
| F: Hermitian dilation | (i) wrong physics (folded picture, linear-T entropy); (ii) same 1/g via σ-gap | — | skip (keep as thesis argument) |

---

## 7. Minimal experiments (precise specs, none run)

- **E1 — phase-rigidity scan (do first, zero theory risk).** For Alcaraz p=0.1 and XXZ Δ=1, over
  the existing T-ladders: after convergence, record rⱼ = |⟨Lⱼ|Rⱼ⟩_noconj| (pre-bi-normalization,
  with ‖Lⱼ‖=‖Rⱼ‖=1) for the cluster members, vs T. Flat O(1) → avoided crossing (proceed as
  planned); →0 → near-EP (switch all cluster handling to the Schur-form variants). Data likely
  regenerable from cached warm-start ladders in a few short runs.
- **E2 — wall-position vs ε sweep (how much of the wall is implementation?).** Alcaraz p=0.1,
  T = 8…12, ΔT=1: rerun the entropy with (i) `cutoff` ∈ {1e-10, 1e-12, 1e-14} in the de-mixing/
  apply schedule, (ii) `trunc_mode=:rtm` vs `:rdm`, (iii) `basis=:eig` vs `:schur` — 3×2×2 grid,
  warm-started. Observable: the T where the dome peak departs >10% from the CFT chord. If the wall
  T shifts with cutoff or trunc_mode, that fraction of the wall is reducible ε; if immobile, it is
  pure g(T) physics and only A/B can help.
- **E3 — refined-Ritz extraction.** Same grid points as E2 at one fixed schedule: after
  convergence, extract (L₀,R₀) both ways (current `eigen(W)` route vs §2.2 minimizer at the tracked
  θ_phys), compare the two entropies and the two residual norms ‖(E−θ)x‖. Prediction registered in
  §2.2.
- **E4 — candidate A.** (a) *Consistency*: in the clean window (Alcaraz p=0.1, T=6..8), compute
  S₂(T;L) for L ∈ {10², 10³, 10⁴}; must converge to the current single-pair answer as g·L grows
  (g~1e-2 there → already converged at L=10³). (b) *Boundary*: compute the cluster overlaps
  ⟨Lᵢ|b⟩ for b = the structured `fw_tMPS` seed and for 3 random b's; check the physical b gives
  stable weights. (c) *Through the wall*: T=9..12, S₂(T;L) family; check continuity in T at fixed
  L, the Im S₂ plateau vs πc/12, and whether the Re dome stays on the chord; attempt a g(T)·L
  collapse. (d) *XXZ cross-check*: Δ=1, T=3..6 across the T≈4 wall with the 4-member cluster —
  the strongest test since the wall is structural there.
- **E5 — candidate B for XXZ.** Derive/solve the on-site gauge matrix V for the rotated-Néel U(dt)
  (Murg order-1 kernel is the easiest analytically); verify symmetry on the block: ‖Q̃W̃ − W̃Q̃‖ on
  the k×k restriction < 1e-6; check Q̃'s cluster eigenvalues split ±1 cleanly at T=4..6; recompute
  per-sector entropies and compare the intra-sector gap to the k=6 |θ₅/θ₁| result (NB9 §4b) as an
  independent confirmation.
- **E6 — ± partner probe (Alcaraz).** Build the temporal staggering operator D = diag(+,−,+,−,…) on
  the time legs (single-site operator product); compute |⟨L₁|D|R₀⟩| and |⟨L₀|D|R₀⟩| at T=6..10. If
  the first is O(1) and the second ~0, the ± pair is D-resolved → candidate B applies to Alcaraz
  too, and the open "why does the −λ₀ partner exist" question (next_steps.md) is answered as
  temporal band-folding.
- **E7 — Procrustes stress test.** Down the existing warm-started ladder: m×m Procrustes-aligned
  frame, carry the clean-window physical vector, entropy vs T; compare drift against attempt #2's
  recorded 0.26→0.34→0.86.

---

## 8. Open questions for the discussion session

1. **Is the wall's entropy even the right target past g·L_phys ~ 1?** Candidate A reframes the
   thermodynamic-limit entropy as undefined at the degeneracy — should the thesis *claim* that
   (with §5.4's argument), making S₂(T;L) the honest observable, rather than presenting any
   technique as recovering "the" T→∞ entropy?
2. **What is the correct physical `b`** for a power-method calculation that never had a finite
   width — is the structured `fw_tMPS` seed the right object, and does its known Z₂-sector
   structure match the weights a finite-chain TDVP calculation would imply?
3. **Does the MPO symmetry-gauge V exist cleanly for the VD2 kernel** (whose enlarged 1+χ+χ²
   virtual leg is the temporal physical leg), or only for the Murg construction? If VD2's V is
   ugly, candidate B runs on the Murg tMPO (construction-independence is already established,
   NB12 §4e).
4. If E1 finds near-EP behavior, the ± pair merging is a Jordan-block formation — is *that* the
   emergent-dual-unitarity endpoint (a physical statement), and is the Puiseux √ scaling visible
   in the θ-splitting data already cached?
5. Cheap-but-unexplored: does running the *final few* iterations (not the whole run) at
   `cutoff=1e-14`/larger χ — after warm-start convergence at the production schedule — reduce the
   effective ε enough to matter, given the wall's ε/g ~ 1 balance? (A one-parameter variant of E2.)

---

## Web sources consulted

- [Refined CJ–SS–RR method / refined Ritz vectors for clustered values (arXiv:2605.12846)](https://arxiv.org/pdf/2605.12846)
- [Jia & Stewart, convergence of Ritz values, Ritz vectors, and refined Ritz vectors](https://www.researchgate.net/publication/2450993_On_the_Convergence_of_Ritz_Values_Ritz_Vectors_and_Refined_Ritz_Vectors)
- [Jia, residuals of refined projection methods (Comput. Math. Appl.)](https://www.sciencedirect.com/science/article/pii/S0898122100003217)
- [Stewart, analysis of Rayleigh–Ritz for eigenspaces](https://www.researchgate.net/publication/2424724_An_Analysis_of_the_Rayleigh--Ritz_Method_for_Approximating_Eigenspaces)
- [Harmonic Ritz convergence analysis (arXiv:1603.01785)](https://arxiv.org/pdf/1603.01785)
- [Homotopy continuation of non-Hermitian dispersion curves through EPs (arXiv:2605.15089)](https://arxiv.org/pdf/2605.15089)
- [Continuous decompositions and coalescing eigenvalues for parameter-dependent matrices](https://www.researchgate.net/publication/282766494_Continuous_Decompositions_and_Coalescing_Eigenvalues_for_Matrices_Depending_on_Parameters)
- [Exceptional points in non-Hermitian systems (overview)](https://www.emergentmind.com/topics/exceptional-points-eps)
- [Perturbation theory in the complex plane: exceptional points (IOP)](https://iopscience.iop.org/article/10.1088/1361-648X/abe795)
- [Zauner-Stauber et al., transfer matrices and excitations with MPS (NJP 17, 053002)](https://iopscience.iop.org/article/10.1088/1367-2630/17/5/053002)
- [Z_N symmetry breaking in PEPS (arXiv:1703.04137)](https://arxiv.org/pdf/1703.04137)
- [Long-range order and symmetry breaking in PEPS (arXiv:1505.04217)](https://arxiv.org/pdf/1505.04217)
- [Infinite boundary conditions / limit cycles in iDMRG (arXiv:1804.09163)](https://arxiv.org/pdf/1804.09163)
