# 5.4 A block power method for the leading transfer spectrum

*(drop-in section for the thesis; written to follow §5.1–5.3, whose notation — the spatial
transfer matrix $E$, its biorthogonal eigenpairs $\langle L_i|$, $|R_i\rangle$ with eigenvalues
$\mu_i$ ordered by modulus — it reuses. Equation numbers are local, (B1)–(B9).)*

## Why a single dominant vector is not enough

The power method of Section 3.3 converges to the dominant eigenpair at a rate set by the
gap ratio $|\mu_1/\mu_0|$: each application of $E$ suppresses the subleading components of
the boundary states by one factor of $|\mu_1/\mu_0| < 1$. Section 5.3 showed, however, that a
quench to a critical point drives exactly this ratio to one — emergent dual unitarity *is*
the statement that the leading moduli become degenerate as $T$ grows. The better the
conformal physics, the worse the power method: at the temporal extents where the
interesting scaling sets in, the single-vector iteration slows down critically and eventually
stalls, mixing the dominant eigenvector with its near-degenerate partners.

The standard cure, familiar from Hermitian eigensolvers, is to iterate a *block* of $k$
vectors at once and let them span the whole quasi-degenerate cluster. Convergence of the
leading pair is then governed not by $|\mu_1/\mu_0|$ but by $|\mu_k/\mu_j|$ — the gap
*between the block and the rest of the spectrum* — which stays finite as long as the block
is large enough to contain the cluster. The subtlety is that our $E$ is non-Hermitian and
non-normal (Section 5.1), so every ingredient of the block method — orthogonalization,
projection, even the notion of "overlap" — must be reformulated in the biorthogonal,
*bilinear* language. This section walks through that reformulation as implemented in
`block_transfer_eigs` (our addition on top of ITransverse.jl).

## The two-sided block and the oblique projection

We keep $k$ right vectors $|R_1\rangle,\dots,|R_k\rangle$ and, independently, $k$ left vectors
$\langle L_1|,\dots,\langle L_k|$, all of them temporal MPS on the same time sites, seeded with
random complex tensors (Section 3.3 explains why a random seed is required). One iteration
applies the transfer matrix to every member of the right block and its *transpose* to every
member of the left block,

$$|\tilde R_j\rangle = E\,|R_j\rangle, \qquad
  \langle \tilde L_j| = \langle L_j|\,E
  \quad\Longleftrightarrow\quad
  |\tilde L_j\rangle = E^{\mathsf T}|L_j\rangle .
  \tag{B1}$$

The transpose — implemented as a prime-level swap of the MPO's physical legs, with *no*
complex conjugation — is the adjoint with respect to the bilinear form of Section 5.1: it is
defined by $(E^{\mathsf T}L)^{\mathsf T}R = L^{\mathsf T}(ER)$, so left overlaps computed with
`overlap_noconj` remain consistent. Using the Hermitian adjoint $E^\dagger$ here would be a
conceptual error, not a convention choice: it would converge the left block to the dominant
eigenvector of $E^\dagger$, which for a non-normal $E$ is *not* the dominant left eigenvector
of $E$.

As in the single-vector method, each MPO application multiplies the bond dimension, so
every product in (B1) is truncated back to a cap $\chi_{\max}$ (more on the truncation mode
below), and each vector is renormalized.

## The projected pencil: Rayleigh–Ritz without conjugation

After the applications, we compress all the spectral information the blocks carry into two
small $k\times k$ matrices,

$$S_{ij} = \langle L_i | R_j\rangle, \qquad
  M_{ij} = \langle L_i | E | R_j \rangle ,
  \tag{B2}$$

both evaluated with the bilinear overlap (no conjugation anywhere — $S$ is complex
*symmetric-free*: it has no Hermiticity properties at all, and that is fine). $S$ is the metric
of the oblique (Petrov–Galerkin) projection: it says how the left test space and the right
trial space are angled against each other. $M$ is the transfer matrix seen through that
projection. Approximations $\theta$ to the eigenvalues of $E$ — *Ritz values* — then solve the
small generalized eigenproblem

$$M\,v = \theta\, S\, v .
  \tag{B3}$$

Numerically we do **not** hand (B3) to a generalized eigensolver. Near the degeneracy the
columns of the right block begin to coincide, $S$ becomes nearly singular, and
`eigen(M, S)` responds with infinite or wildly scattered eigenvalues. Instead we regularize
with the pseudoinverse,

$$W = S^{+} M, \qquad W = V\,\Theta\,V^{-1},
  \tag{B4}$$

where $S^{+}$ is the Moore–Penrose pseudoinverse with a relative tolerance ($10^{-12}$):
directions of the block that have effectively collapsed onto each other are projected out
rather than inverted. The eigenvalues $\Theta = \mathrm{diag}(\theta_1,\dots,\theta_k)$, sorted
by modulus, are our estimates of $\mu_0,\dots,\mu_{k-1}$; the columns of $V$ say which linear
combination of the current right block best approximates each eigenvector.

## Pairing the left coefficients exactly

The left block needs its own mixing coefficients, and here hides a latent trap. One can
diagonalize the transposed pencil, $u^{\mathsf T} M = \theta\, u^{\mathsf T} S$, in a second,
independent `eigen` call — but then nothing guarantees that its eigenvalues come out in the
same order, and inside a quasi-degenerate cluster "match each left eigenvalue to the nearest
right one" is ambiguous precisely when it matters: several $\theta$'s sit within a fraction
of a percent of each other and the greedy match *can* pick a wrong partner. (An honest
aside: controlled dense tests with planted clusters showed that in exact arithmetic this
mispairing does **not** in fact destabilize the iteration — the subspace self-corrects — and
the noisy tricritical gap data we once blamed on it turned out to be something else entirely:
a *selection* artifact, the "physical $\lambda_0$" rule hopping between members of a genuinely
near-degenerate band. The lesson is sharper than the bug hunt: inside a band, do not select a
single eigenvalue — report the band. The exact pairing below is still the right construction:
it is cheaper, and it eliminates the latent failure mode by design.)

Indeed, no second diagonalization is needed. Writing $y^{\mathsf T} = u^{\mathsf T} S$,
the left eigenproblem becomes $y^{\mathsf T} W = \theta\, y^{\mathsf T}$: the $y$'s are the left
eigenvectors of the *same* matrix $W$, i.e. the rows of $V^{-1}$. Hence

$$u_j = \big(S^{+}\big)^{\mathsf T} \big(V^{-1}\big)^{\mathsf T} e_j ,
  \tag{B5}$$

paired with $\theta_j$ *by construction*, and automatically biorthogonal through the metric,
$u_i^{\mathsf T} S\, v_j = \delta_{ij}$ (up to the pseudoinverse regularization). One
decomposition, exact pairing, no heuristic.

## De-mixing, and why the eigenvector basis eventually fails

The iteration closes by rotating the applied blocks onto the Ritz combinations,

$$|R_j\rangle \leftarrow \sum_i V_{ij}\, |\tilde R_i\rangle, \qquad
  |L_j\rangle \leftarrow \sum_i U_{ij}\, |\tilde L_i\rangle ,
  \tag{B6}$$

(the MPS linear combinations are formed with an exact direct-sum addition and compressed
once at the end — the naive density-matrix addition is numerically fragile for
near-parallel vectors). Each new pair is then truncated, in one of two modes that mirror
the RTM/RDM discussion of Section 3.2:

- **RTM mode** truncates the matched pair $(L_j, R_j)$ *jointly* on its bilinear transition
  matrix $|R_j\rangle\langle L_j|$ — optimal for the physical overlap, dramatically cheaper in
  bond dimension, but built on the SVD of a non-Hermitian object that becomes
  ill-conditioned exactly at the degeneracy;
- **RDM mode** truncates every vector *independently* on its Hermitian density matrix
  $|v\rangle\langle v^*|$ — blind to the $L$–$R$ pairing, more expensive, but positive and
  well-conditioned through the cluster. It is the safe fallback when RTM output turns noisy.

There remains one structural weakness that no truncation mode removes. Equation (B6)
de-mixes onto the *eigenvector* basis, and for a non-normal matrix the eigenvector basis
itself degenerates as eigenvalues cluster: the condition number of $V$ grows like the
inverse gap. Rotating the block onto nearly-parallel directions amplifies every source of
noise (truncation, finite convergence) by that condition number. The remedy, standard in
non-Hermitian subspace iteration, is to give up on individual eigenvectors *inside* the
cluster and maintain only a well-conditioned basis of the invariant *subspace*: we
QR-orthonormalize the modulus-sorted coefficient matrices,

$$V \to Q_V, \qquad U \to Q_U, \qquad (Q^\dagger Q = \mathbb 1),
  \tag{B7}$$

an ordered-Schur-like choice whose leading $j$ columns span the same leading-$j$
eigenspaces while the rotations stay unitary. The Ritz values are still read from the pencil
(B3)–(B4) at every step; only the basis carried between iterations changes. In this mode the
individual pairs $(L_j, R_j)$ for $j$ inside a cluster are no longer eigenvector
approximations — a price one must remember when feeding them into entropy formulas —
but the *spectrum*, and the leading pair when it is isolated, remain clean.

## Convergence, the $\pm$ partner, and the physical gap

Convergence is declared on the movement of the leading Ritz values between iterations,
$\max_j |\theta_j^{(n)}-\theta_j^{(n-1)}| < \varepsilon$. One more non-Hermitian subtlety
enters here: the transfer spectrum comes in near-$\pm$ pairs — for each dominant $\mu_0$
there is a partner of nearly equal modulus and phase shifted by $\approx\pi$ (a
period-two structure in the time direction). Two consequences:

1. *Bookkeeping.* Sorting by modulus lets $\mu_0$ and its partner swap places between
   iterations, which a naive difference registers as a spurious jump; the tracked values must
   be matched between iterations by continuity in the complex plane before differencing.
2. *Physics.* The naive gap ratio $|\theta_2|/|\theta_1|$ measures the *splitting of the pair*,
   not the spectral gap. The meaningful, "physical" gap takes the largest modulus among the
   Ritz values that are **neither** the physical $\mu_0$ (selected by continuity in $T$)
   **nor** its $-\mu_0$ partner:

$$\mathrm{gap}(T) = \frac{\max\{|\theta_i| : i \neq i_0,\, i \neq i_{\text{partner}}\}}{|\theta_{i_0}|}.
  \tag{B8}$$

   All gap-closing plots in this thesis use (B8); the RTM- and RDM-truncated results agree on
   it, which excludes a truncation artifact.

Two practical accelerations complete the algorithm. Along a ladder of temporal extents
$T_1 < T_2 < \dots$ the converged block at one rung *warm-starts* the next: the converged
tensors are re-indexed onto the longer time-site set and the few added sites are padded
with a small random tail, so the iteration only has to relax the new sites. And the bond-cap
and cutoff follow a schedule — loose and cheap in the early iterations, tight near
convergence — since the early block is random garbage not worth representing accurately.

## What bounds the reach

The block method extends the usable window considerably — it is what carried the Alcaraz
analysis through the gap closing to $T \approx 9$ — but it does not abolish the wall of
Section [limits]. The reason is worth stating precisely, because it separates what an
algorithm can and cannot fix. The leading eigen*values* of $E$ are well-conditioned
essentially always (their sensitivity involves the biorthogonal pair through
$1/|u^{\mathsf T} v|$, which the block method keeps under control). The individual
eigen*vectors*, however, have condition number $\sim 1/\text{gap}$: as emergent dual
unitarity closes the gap, any numerical representation of "the" dominant eigenvector —
and hence any entropy built from it — degrades at a rate fixed by the physics itself, not by
the iteration. A block method with exact pairing and a Schur basis therefore delivers
clean eigenvalues (moduli, phases, the gap (B8)) arbitrarily deep into the dual-unitary
regime, while the entropy profile, which needs the eigenvector, is trustworthy only in the
pre-wall window where the leading pair is still isolated:

$$\underbrace{\text{spectrum of } E}_{\text{robust: eigenvalue conditioning}}
  \qquad\text{vs}\qquad
  \underbrace{S^{\mathrm{gen}}_2 \text{ from } (L_0, R_0)}_{\text{bounded by } 1/\text{gap}} .
  \tag{B9}$$

Characterizing this boundary — where it sits for each model, and why frustration and a
larger central charge pull it to earlier times — is part of the results of Sections
[gap-closing] and [tricritical].

One escape from the asymmetric bound deserves mention: when the model admits a left–right
**symmetric** MPO (as Ising does via the Murg construction), the Autonne–Takagi route replaces
the two-sided iteration entirely, with qualitatively better conditioning near degeneracies. For
the rotated XXZ chain such a construction exists — its two-site term decomposes into mutually
commuting $XX$, $YY$, $ZZ$ layers, each an exact bond-2 symmetric factor — and notebook 12 uses
it to measure directly how much of the asymmetric method's early wall is the asymmetry itself.
