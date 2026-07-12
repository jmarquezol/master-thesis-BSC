# The symmetric-Takagi route for XXZ-Néel: the gauge question resolved
### Research report — which symmetry the transverse machinery actually needs, why the σʸ obstruction does not block it, and the direct dial-(iii) verdict
*(companion working notebook: `NBs/12_xxz_symmetric_mpo.ipynb`; closes the "one honest loose end"
of CLAUDE.md §18 and the final caveat paragraph of the barrier section. Written 2026-07-12.)*

---

## 0. Executive summary

The §18 loose end was: *"we could not run a correctly-gauged symmetric-Takagi solver (no
single-site symmetric gauge exists for XXZ)"*. This report shows that statement rests on a
**misidentification of which index swap the symmetric machinery requires**, and that once the two
swaps are disentangled:

1. **The XXZ-Néel Murg tMPO (both orders) is already correctly gauged.** The tMPO built from
   `expH_xxz_neel_murg` is an *exactly* transpose-symmetric operator on its temporal-physical
   legs — dense contraction gives ‖M−Mᵀ‖/‖M‖ = 1.1×10⁻¹⁸ (order 1) and 9.2×10⁻²¹ (order 2) — and
   its left/right transfer eigenvectors coincide to |⟨l,r⟩| = 1.0000 (dense eig, all four leading
   members). The ⟨L| = |R⟩ᵀ identification that `powermethod_sym` + Autonne–Takagi rely on is
   **exactly valid** for this construction.
2. **The σʸ obstruction (NB12 §3b) is real but constrains the *other* index pair** — the one that
   becomes the tMPO's internal *time-bond* direction, which no part of the symmetric pipeline
   (`powermethod_sym`, `RTMsym` truncation, `generalized_vn_entropy_symmetric`) ever transposes.
   The §3b conclusion "the symmetric-Takagi machinery is structurally unavailable to XXZ" is
   **withdrawn**; what remains true is the narrower statement that no single-site gauge makes the
   tensor symmetric under *both* swaps simultaneously (irrelevant to the echo pipeline).
3. **Consequently NB12's §4/§4b/§4c `powermethod_sym` runs are resurrected as the valid, direct
   dial-(iii) experiment.** Their content: the symmetric n→1 dome inflates at the **same T≈4 wall**
   as the asymmetric Rényi-2 dome (Δ=0.5 peaks 0.55→0.86→1.10 over T=3→5; Δ=1.0 jump 0.43→1.03 at
   T=4→5, matching VD2's 0.25→0.81), seed-independent to all printed digits, warm/cold-consistent.
   **The wall does not move. Dial (ii) dominates dial (iii) — now shown directly, not only via the
   gauge-free surrogates (§4d/§4e).**
4. **The package's `expH_XXZ_svd` SymSVD failure is diagnosed**: it sweeps *non-commuting* bond
   gates sequentially, and spatial reflection reverses that order — the asymmetry scales exactly as
   O(dt²) (measured 1.22×10⁻² → 3.06×10⁻³ → 7.65×10⁻⁴ under dt = 0.2 → 0.1 → 0.05; ratios 3.99,
   4.00). It is an assembly-order artifact, not a basis/frame obstruction; the Potts template works
   because Potts's bond gates are diagonal (all commute). The commuting-*layer* factorization
   (`exp2site_murg`) is the correct XXZ analog of the Potts trick, and it already exists.
5. A new observable closes the conditioning question: **phase rigidity r_j(T)** of the leading
   transfer eigenvectors, symmetric-Murg vs asymmetric-VD2 vs Ising (cache
   `results/data/nb12_rigidity.jld2`, §6 + UPDATE below). Answer: the symmetric construction
   **relabels** the near-exceptional-point collapse, it does not avoid it — sym and asym
   rigidities fall at the same geometric rate (constant factor ≈2–3, same exponent), with the
   symmetric arm's L/R alignment pinned at 1.000 (E=Eᵀ live-confirmed). And a genuine surprise
   from the Ising control: **its rigidity collapses even deeper (r₀→1.3×10⁻³ by T=12) while its
   entropy stays clean to T≈14** — so r→0 is the universal dual-unitarity fingerprint, *not* the
   wall criterion; the wall is set by *when* the modulus band tightens to few-% (XXZ: T≈4, driven
   by its exactly degenerate Z₂ pair; Ising: T≳12, smooth non-degenerate band).

**Verdict (item 4 of the task):** symmetry does *not* fix the XXZ block-PM wall. The correctly
gauged symmetric-Takagi solver was, in fact, already run — it walls at T≈4 exactly like the
two-sided solver, because the exact Néel Z₂ degeneracy is a property of the quench's transfer
matrix (present in the symmetric tMPO's own spectrum, §4d) and eigenvector conditioning inside an
exactly degenerate cluster is basis-independent-ly ill-posed. The §18 loose end is closed in the
expected direction, but now by a *valid direct experiment* rather than by expectation.

---

## 1. The two symmetries of a rotated tMPO tensor (the disentanglement)

The 90° rotation maps the spatial U(δt) MPO tensor W[α, β; s, s′] (α,β = spatial links, s,s′ =
spatial physical) onto the tMPO tensor as follows (verified against the installed
`make_fwtmpoblocks`, `~/.julia/packages/ITransverse/8pmYI/src/tmpo/fw_tmpo_blocks.jl:90-101`):

| spatial legs of W | become | tMPO role | checker line | XXZ-Murg status |
|---|---|---|---|---|
| links α, β (`iLink1`,`iLink2`) | → | **temporal-physical** `Site,time` pair (iP, iPs) — the legs the tMPS attaches to | "bond(space) => phys(time)" | **symmetric, 4×10⁻¹⁹** |
| physical s, s′ (`icP`,`icP′`) | → | **temporal links** (iL, iR) — internal bonds of the tMPO chain | "physical(space) => bond(time)" | *fails, 0.45 (σʸ)* |

The temporal transfer operator E (the tMPO as a matrix acting on boundary tMPS vectors) has matrix
elements E[{α_t},{β_t}] = ∏_t W[α_t, β_t; ·] with the s-legs contracted internally along the time
direction. Its transpose swaps {α_t}↔{β_t} slice by slice, so

> **E = Eᵀ ⟺ the bulk tensor is invariant under the spatial-link swap α↔β** (plus the same for
> the imaginary-time blocks, and time-boundary attachments that do not touch the α,β legs — both
> verified: `fw_tMPO` contracts `bl`/`tr` on the outer *time-link* legs only,
> `build_fw_tmpo.jl:221-239`).

The spatial-*physical* swap s↔s′ is untouched by this transpose: it corresponds to
**time-reflection** of the column (equivalently, to operator transpose symmetry U=Uᵀ of the
spatial propagator, up to link gauge). Geometrically: transposing E mirrors the 2D echo network
across a *vertical* axis (space reflection); the σʸ obstruction lives on the *horizontal* mirror
(time), which the Loschmidt-echo fixed-point equations never perform.

Three operator-level properties are therefore mutually independent, and the numerical experiment
(§3) exhibits all combinations:

| construction | reflection PUP=U (operator) | transpose U=Uᵀ (operator) | **E=Eᵀ raw (what the solver needs)** |
|---|---|---|---|
| Ising Murg | ✓ | ✓ | ✓ (0.0) |
| XXZ Murg, order 2 (palindrome) | ✓ | ✓ (3×10⁻¹⁸) | **✓ (9×10⁻²¹)** |
| XXZ Murg, order 1 | ✓ | ✗ (1×10⁻³, not a palindrome) | **✓ (1×10⁻¹⁸)** |
| XXZ VD2 | ✓ | ✓ (7×10⁻²⁰) | ✗ (0.35) |
| package `expH_XXZ_svd` staircase | ✗ O(dt²) | ✗ O(dt²) | ✗ |

Note in particular the order-1 row: an operator-*asymmetric* propagator whose tMPO is *exactly*
symmetric — and the VD2 row: an operator-symmetric propagator whose tMPO is not. "Symmetric model"
was never about the operator; but it was also never about the *time-direction* tensor gauge — it is
about the **spatial-reflection gauge being trivial (G = 1) in the raw index basis**. The subtle
point that makes "manifest" non-negotiable: for an injective reflection-symmetric MPO, the
reflected tensor always satisfies W^refl = G W G⁻¹ for some link gauge G (fundamental theorem);
but after rotation the spatial links are *promoted to the temporal physical basis*, so G is no
longer an internal gauge — it is a similarity transformation on the temporal Hilbert space. E and
Eᵀ are then merely similar, not equal, and a raw-basis `powermethod_sym` is wrong by exactly that
G-twist (VD2's case). A manifestly symmetric construction is one with G = 1; Murg-type layer
factorizations deliver G = 1 automatically because both halves of every layer factor are the
*same* operator (`exp2site_murg`: W[1,2] = W[2,1] = √(i sin)√(cos)·σᵃ).

### 1.1 What the installed pipeline actually assumes (code audit)

Every stage of the symmetric chain was read in the installed source
(`~/.julia/packages/ITransverse/8pmYI/src/`):

- `powermethod_sym` (`power_method/symm_pm.jl`): never constructs ⟨L|; applies the tMPO via
  `tapply` and normalizes by `overlap_noconj(psi, psi)`. Sole assumption: the left fixed point is
  the unconjugated transpose of the right one ⟺ E = Eᵀ on the `Site,time` legs.
- `RTMsym` truncation (`truncation_sweeps/sweeps_sym.jl`): all environments are ψ doubled against
  its own *primed, unconjugated* copy (`env *= Ai; env *= noprime(Ai′,…)`; `rho = E[j+1]*L*L″`).
  No `dag` anywhere; no time-bond transposition.
- `generalized_vn_entropy_symmetric` / `diagonalize_rtm_symmetric` (`entropies/`): same bilinear
  doubling (`env *= psi[ii]; env *= psiP[ii]`), normalization `overlap_noconj(psi,psi)`.

Nothing in the chain ever swaps the tMPO's link legs (the σʸ-obstructed pair). The
`FwtMPOBlocks` checker tests both swaps and *warns* on either; only the "bond(space) => phys(time)"
leg is load-bearing for this pipeline. (The "physical(space) => bond(time)" leg matters for
*folded* / time-reflection constructions — `fold_tmpo_blocks.jl` — which the echo pipeline does not
use.)

### 1.2 The exact condition (task item 1) and where the σʸ argument survives

**Condition for the symmetric contraction of a uniform tMPO:** the spatial U(δt) MPO's bulk tensor
must be invariant under the raw spatial-link swap. Given operator reflection symmetry, this is a
gauge condition (G = 1 above); a manifestly symmetric form exists iff the reflection gauge G can be
chosen trivial — equivalently (when G is symmetric) absorbed by conjugating the links with G^{1/2}.

- *Single same-operator layers* exp(θ Σ σᵃσᵃ): Murg split gives equal halves ⟹ G = 1. Works for
  **any** Pauli, σʸ included — transpose-antisymmetry of σʸ is invisible to the link swap because
  both link entries store the *same* operator.
- *Products of such layers*: the link swap acts layer-wise on the stacked links and never reorders
  the physical operator string, so any product of individually link-symmetric layers is
  link-symmetric — **regardless of layer order** (this is why order 1 works despite U ≠ Uᵀ). The
  palindrome buys Trotter order, not symmetry (as NB12 §1b already noted).
- The **σʸ obstruction** (NB12 §3b) is correct mathematics aimed at the *other* pair: no single-site
  basis makes σˣ, σʸ, σᶻ simultaneously transpose-symmetric, so the XXZ tensor can never pass
  *both* checker legs at once. Standing consequence: the *folded*/time-reflected machinery (and
  any future use of the time-bond transpose) genuinely is unavailable to XXZ in a single-site
  gauge. For the echo pipeline: moot.

---

## 2. Why the package's XXZ SymSVD fails, and why Potts works (task item 1, diagnosis)

`expH_XXZ_svd` (`chain_models/xxzmodel.jl`) copies the Potts template
(`expH_potts_symmetric_svd`): split each two-site gate by `symm_svd` (A = U·s·Uᵀ, halves =
transposes — the primitive is sound, and the XXZ gate qualifies: exp(ε(XX+YY+ΔZZ)) is site-swap
symmetric, so its ((s₁s₁′),(s₂s₂′)) matrix is complex symmetric), then assemble by a **single
left-to-right sweep of overlapping gates** g₁₂·g₂₃·…·g_{N−1,N}.

- For **Potts** (and Ising ZZ/XX): the coupling operators are diagonal/identical, all bond gates
  **commute**, the sweep equals exp(εΣh) exactly, and the assembly is reflection-invariant. Both
  checker legs pass (the halves are diagonal ⟹ also transpose-symmetric).
- For **XXZ**: adjacent bond gates do **not** commute. The sweep is (i) not exp(εΣh) — it is a
  staircase Trotterization — and (ii) not reflection symmetric, because reflection reverses the
  product order: P(∏ₙgₙ)P = ∏ₙgₙ in *descending* order ≠ ascending. Measured (N=4, Δ=0.5): dense
  reflection asymmetry 1.22×10⁻² / 3.06×10⁻³ / 7.65×10⁻⁴ at dt = 0.2 / 0.1 / 0.05 — ratio 4.0 per
  halving, i.e. **exactly the O(dt²·[h,h′]) commutator scale**. The historical "normdiff
  0.07–0.45" (CLAUDE.md §3b) is the same effect read at tensor level (bulk link-swap normdiff
  0.048–0.18 over the same dt range, plus the O(1) transpose-leg failure).

So: the *primitive* (symmetric SVD split of a site-swap-symmetric gate) transfers to XXZ; the
*assembly* does not. The correct XXZ assembly is the **commuting-layer factorization** — exactly
`exp2site_murg`'s XX/YY/ZZ layers, each internally commuting so each layer is exact and
reflection-manifest, sandwiched into U(δt). That construction exists, is verified against TDVP
(NB12 §1–2), and — per §1 — is correctly gauged for the symmetric solver. An even/odd-brick
symm_svd Trotterization would also be reflection-symmetrizable at the operator level but breaks
the single-site uniformity `make_fwtmpoblocks` requires (`@assert length(eH) == 3`), so the layer
route is also the only one compatible with the package's tMPO builder.

---

## 3. Numerical falsification of the gauge claim (task item 2 — no new construction needed)

Script: session scratchpad `gauge_check.jl` (results reproduced as a cell in NB12 §3c; everything
below is exact/dense, no power method). Predictions P1–P6 all confirmed, Δ=0.5, dt=0.05, nbeta=0:

```
[tensor] Ising-Murg      swap(iP,iPs)=0.0        swap(iL,iR)=5.3e-17
[tensor] XXZ-Murg(1)     swap(iP,iPs)=4.4e-19    swap(iL,iR)=3.1e-01   <- σy, harmless
[tensor] XXZ-Murg(2)     swap(iP,iPs)=2.5e-19    swap(iL,iR)=3.1e-01
[tensor] XXZ-VD2         swap(iP,iPs)=3.8e-01    swap(iL,iR)=3.4e-01   <- genuinely asymmetric
[dense]  Ising-Murg      Nt=3 d_t=2    |M-M^T|/|M| = 0.0
[dense]  XXZ-Murg(1)     Nt=3 d_t=8    |M-M^T|/|M| = 1.1e-18          <- E = E^T exact
[dense]  XXZ-Murg(2)     Nt=2 d_t=32   |M-M^T|/|M| = 9.2e-21
[dense]  XXZ-VD2         Nt=3 d_t=13   |M-M^T|/|M| = 3.5e-01          <- control
[eig]    XXZ-Murg(1): all 4 leading members  |<l,r>| = 1.0000, r_j = selforth (complex-symmetric)
[eig]    XXZ-VD2:     |<l,r>| = 0.969/0.989/0.989/0.935 (left ≠ right)
```

The imaginary-time (nbeta) blocks pass the same link-swap check (4×10⁻¹⁹-level; also visible in
NB12's §4 `FwtMPOBlocks` logs, which show every "bond(space) => phys(time)" line as *Info: Tensor
symmetric* while only the "physical(space)" lines warn). At the tiny dense-checkable T the phase
rigidities of Murg and VD2 are nearly identical (r₀ = 0.977 vs 0.969) — the constructions start
from the same conditioning; what matters is their behavior at the wall (§6).

---

## 4. The resurrected direct experiment (task item 3, part 1)

With the gauge objection removed, NB12 §4/§4b/§4c are re-classified from "cautionary record" to
**the valid direct dial-(iii) test** — `powermethod_sym` (alg `RTMsym`, opt `:RTM_R`) + Takagi n→1
entropy on the exactly-symmetric Murg tMPO. What that data says (caches `nb12_xxz_sym.jld2`,
`nb12_sym_seedtest.jld2`, `nb12_xxz_sym2.jld2`):

| Δ | T | sym χ | sym n→1 peak | asym VD2 S₂ peak (NB9) |
|---|---|---|---|---|
| 0.5 | 2 | 5 | 0.5625 | 0.2179 |
| 0.5 | 3 | 8 | 0.5472 | 0.3891 |
| 0.5 | 4 | 8 | **0.8632** | **0.7023** |
| 0.5 | 5 | 10 | **1.0963** | 0.8248 |
| 0.5 | 6 | 15 | 1.1478 | 0.9723 |
| 1.0 | 4 | 10 | 0.4315 | 0.2500 |
| 1.0 | 5 | 13 | **1.0297** | **0.8094** |
| 1.0 | 6 | 17 | 1.0131 | 0.8891 |

(Absolute values differ between columns by construction — n→1 with coefficient c/6 vs Rényi-2 with
c/8 and different s₀ — the comparison is the *location of the inflation*.) The symmetric dome
inflates between T=3→5 (Δ=0.5) and T=4→5 (Δ=1.0) — **the same wall, at the same T, as the
asymmetric route**; past it the dome sits >1 and chord fits return noise (sym "c" column of the
NB12 §4c table swings 0.08–1.48 with no clean window — consistent with a wall at T≈4 leaving only
T≤3, too short for a slope). The §4b seed test is exact to all printed digits across three cold
random seeds at T=4 and T=6 — like the asymmetric route, the truncated symmetric iteration has a
unique deterministic attractor past the wall; and §4c's warm-started ladder reproduces phase 1's
cold-started points, eliminating the cold-start hypothesis. The Takagi route's failure signature is
its own (RTM "norm² not real" warnings for T≳5 — quasi-null Takagi vectors), but its *location* is
unchanged.

Two independent confirmations already in NB12 stand unchanged and are now corroborating rather
than load-bearing: the spectrum bridge (§4d — same 4-fold Z₂ band in the symmetric tMPO at the
same T≈4) and the construction-independence test (§4e — Murg through the two-sided solver
reproduces the VD2 dome to 3–4 digits including the inflation).

---

## 5. What this does to the barrier narrative

- **"Dial (ii) dominates dial (iii)" is now a direct experimental statement.** The strongest
  version of dial (iii) — a manifestly symmetric MPO contracted by `powermethod_sym` with the
  Autonne–Takagi RTM, the exact machinery that carries Ising to T=14 — was validly applied to XXZ
  and does not move the T≈4 wall by any amount.
- **The "structurally unavailable" claim is withdrawn** from the thesis narrative (barrier section
  §4's parenthetical loose end): the symmetric-Takagi route *is* available to XXZ. Its
  unavailability was the one thing keeping the dial-(iii) conclusion indirect. What survives of the
  σʸ analysis: XXZ's tensor cannot be symmetric under both swaps at once, so Ising remains
  *doubly* privileged (its X-coupling is transpose-symmetric, making even the time-bond leg
  manifest) — but that extra privilege is irrelevant to the echo contraction.
- **Ising's T=14 reach is re-attributed**: not "enabled by having a symmetric MPO at all" (XXZ has
  one too, and it buys nothing at the wall) but by *what the symmetric machinery is conditioning
  against* — an emergent, slowly-tightening near-degeneracy with no exact symmetry cluster on top,
  plus small d_t and c=1/2. Symmetry improves constants where the eigenvector question is
  well-posed; it cannot make an ill-posed question well-posed. (Formal version: every square
  complex matrix is similar to a complex-symmetric one, so complex symmetry per se imposes *no*
  constraint on eigenvector conditioning; near-EP self-orthogonality is fully compatible with
  M = Mᵀ, where it appears as Takagi quasi-null vectors, vᵀv → 0.)

## 6. Phase rigidity across the wall: does symmetry delay the collapse?

The wall mechanism identified in `eigvec_robustness_report.md` is a near-exceptional-point: the
bilinear condition number κ(λ₀) = ‖L₀‖‖R₀‖/|⟨L₀|R₀⟩| grows geometrically (Alcaraz p=0.1: 12→1413
over T=2..6, i.e. phase rigidity r₀ = 1/κ collapsing 0.085→0.0007). The sharp remaining question:
does the *symmetric* XXZ construction — where L_j ≡ R_j and r_j becomes the Takagi
self-orthogonality |v_jᵀv_j|/‖v_j‖² — avoid or merely relabel that collapse?

Measurement: `block_transfer_eigs` (the gauge-free two-sided solver, exact-diagonalization-validated
in NB5) on three arms, per-member r_j = |⟨L_j|R_j⟩|/(‖L_j‖‖R_j‖) plus the Hermitian alignment
|⟨L_j,R_j⟩|/(‖L_j‖‖R_j‖) (= 1 iff L_j = R_j up to phase, a per-T live check of E=Eᵀ):

- `(:xxzsym, T)` — Murg(1) tMPO, Δ=0.5, dt=0.05, nbeta=4 (the §4d configuration), T=1..6;
- `(:xxzasym, T)` — VD2 tMPO, same physics, T=1..6;
- `(:ising, T)` — Ising Murg control, dt=0.1, nbeta=4, T=2..12.

Cache: `results/data/nb12_rigidity.jld2` (crash-safe, RAM-gated driver; log
`results/logs/nb12_rigidity.log`). **Results: see the table in NB12's rigidity section — filled in
from the cache by its analysis cell.** Interpretation key, pre-registered: if r_j(sym) stays O(1)
through T≈4 while r_j(asym) collapses, symmetry genuinely improves the conditioning and the T≈4
wall must then be blamed on truncation alone; if both collapse together at T≈4, the near-EP is a
property of the quench's transfer matrix in *any* gauge, and the symmetric route's identical wall
(§4) is fully explained.

*(Section finalized from the ladder output — see NB12 and the UPDATE below.)*

---

## UPDATE (2026-07-12, post-ladder): rigidity results — collapse is gauge-blind, and it is NOT the wall criterion

Full ladder (cache `nb12_rigidity.jld2`; k=4, maxdim=48, all points `converged` except Ising
T=10,12 `stuck`, whose values are still indicative):

| T | r₀ sym | r₀ asym | sym/asym | \|θ\| (identical both constructions, 4 digits) |
|---|---|---|---|---|
| 1 | 0.662 | 0.544 | 1.2 | [0.920, 0.322, 0.322, 0.166] |
| 2 | 0.446 | 0.301 | 1.5 | [0.797, 0.588, 0.588, 0.524] |
| 3 | 0.335 | 0.190 | 1.8 | [0.837, 0.769, 0.767, 0.767] |
| 4 | 0.222 | 0.102 | 2.2 | [0.864, **0.829, 0.829**, 0.826] |
| 5 | 0.130 | 0.049 | 2.6 | [0.853, **0.827, 0.827**, 0.827] |
| 6 | 0.092 | 0.029 | 3.2 | [0.856, **0.838, 0.838**, 0.833] |

Ising control (dt=0.1): r₀ = 0.303 / 0.102 / 0.035 / 0.012 / 0.004 / 0.0013 over T = 2 / 4 / 6 /
8 / 10 / 12. Sym-arm Hermitian alignment |⟨L_j,R_j⟩|/(‖L‖‖R‖) = 1.000 at every T for the three
leading members (E=Eᵀ live-confirmed through the neutral solver); asym alignment falls to ~0.01.

**Three conclusions:**

1. **The collapse is gauge-blind.** Sym and asym rigidities fall at the same geometric rate
   (~×0.65–0.7 per unit T); the manifest symmetry buys a slowly growing constant factor (1.2→3.2
   over T=1..6), not a change of exponent. The pre-registered "both collapse together" outcome
   holds: the near-EP is a property of the quench's transfer matrix in any gauge, fully explaining
   why the valid symmetric-Takagi run (§4) walls at the same T.
2. **Rigidity collapse per se is NOT the wall criterion — the Ising control proves it.** Ising's
   r₀ falls *deeper* than XXZ's (1.3×10⁻³ by T=12) while its entropy famously stays on the CFT
   chord to T≈14. Consistent with NB13's Alcaraz observation (r₀ ≈ 0.02 at T=3, dome clean to
   T≈9), the r→0 collapse is the universal emergent-dual-unitarity fingerprint, present in every
   model, symmetric or not.
3. **What actually sets the wall is *when* the modulus band arrives.** XXZ carries an *exactly*
   degenerate pair (|θ₂| = |θ₃| to all printed digits at every T, both constructions — the two
   Néel Z₂ sectors) that closes on λ₀ to ~4% already at T=4, exactly where both domes inflate;
   Ising's four members stay mutually *non-degenerate* and reach comparable few-% tightness only
   at T≳12 — just before its own T≈14 wall. The wall criterion is band tightness (few %), whose
   arrival time is pure quench physics; inside an exactly degenerate subspace the individual
   eigenvector is not ill-conditioned but *undefined*, and no gauge, construction, or solver can
   define it. This is the cleanest quantitative form of "dial (ii) dominates dial (iii)" the
   campaign has produced.

---

## 7. Files touched / reproduction

- `NBs/12_xxz_symmetric_mpo.ipynb` — §3b rewritten (two-swaps disentanglement), new §3c
  (falsification tests), §4/§4b/§4c banners corrected (valid, resurrected), new rigidity section,
  §5 verdict rewritten.
- `results/data/nb12_rigidity.jld2` — new cache (rigidity ladder; regenerable from NB12's cell).
- No library changes: `src/models.jl`'s `exp2site_murg`/`expH_xxz_neel_murg` were already correct;
  the installed ITransverse machinery is used as-is.
- Dense falsification: NB12 §3c cell (originally session scratchpad `gauge_check.jl`).

Key sources verified in the installed package (`~/.julia/packages/ITransverse/8pmYI/src/`):
`tmpo/fw_tmpo_blocks.jl` (rotation + checker), `tmpo/build_fw_tmpo.jl` (boundary attachment),
`power_method/symm_pm.jl`, `truncation_sweeps/sweeps_sym.jl`, `entropies/diagonalize_sym_rtm.jl`,
`entropies/gen_sym_entropies.jl`, `ITenUtils/svd_sym.jl`, `ITenUtils/itensor_utils.jl`,
`chain_models/potts.jl`, `chain_models/xxzmodel.jl`, `chain_models/ising_parallel.jl`.
