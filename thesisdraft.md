| Long-time  |                | dynamics |               | in quantum    |                 | many-body | systems |
| ---------- | -------------- | -------- | ------------- | ------------- | --------------- | --------- | ------- |
| using      | Tensor         |          | Networks      |               |                 |           |         |
| Joaquín    | G. Márquez     | Olguín   |               |               |                 |           |         |
| Supervised | by:            | Stefano  | Carignano     |               |                 |           |         |
| Barcelona  | Supercomputing |          | Center, Plaça | Eusebi Güell, | 08034 Barcelona |           |         |
Aug –, 2026
|     | Brief | description | of the | Master Thesis. |     |     |     |
| --- | ----- | ----------- | ------ | -------------- | --- | --- | --- |
Keywords: TFM, Tensor Networks, Transverse Contraction, Many-body Systems, Quan-
tum
Acknowledgements
I thank my dad, my brother, but specially to my mom, my gf Helena, my flatmates...
| Joaquín | G. Márquez | Olguín: | jmarquol30@alumnes.ub.edu |     |     |     |     |
| ------- | ---------- | ------- | ------------------------- | --- | --- | --- | --- |
1

Contents
| 1 Introduction   | and Motivation |            |             |     | 3   |
| ---------------- | -------------- | ---------- | ----------- | --- | --- |
| 2 The ANNNI-type |                | Model      |             |     | 4   |
| 3 Time Evolution | via            | Transverse | Contraction |     | 6   |
3.1 Temporal entanglement and the CFT connection . . . . . . . . . . . . . . . 7
3.2 Truncation based on RDM and RTM . . . . . . . . . . . . . . . . . . . . . . 8
3.3 The Power Method . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . 9
3.4 Generalized temporal entropies . . . . . . . . . . . . . . . . . . . . . . . . . 10
e−iHt
| 4 MPO Representation |     | of  |     |     | 11  |
| -------------------- | --- | --- | --- | --- | --- |
4.1 The Hamiltonian as a Finite-State Machine . . . . . . . . . . . . . . . . . . 12
4.2 Exponentiating the Finite-State Machine . . . . . . . . . . . . . . . . . . . . 13
5 The Non-Hermitian Transfer Matrix and Dynamical Phase Transitions 15
5.1 Left and right eigenvectors, and biorthogonality . . . . . . . . . . . . . . . . 15
5.2 Eigenvalues versus singular values . . . . . . . . . . . . . . . . . . . . . . . . 15
5.3 Emergent dual unitarity as a spectral statement . . . . . . . . . . . . . . . . 15
| 6 Conclusions  | and Outlook |          |             |            | 16  |
| -------------- | ----------- | -------- | ----------- | ---------- | --- |
| Bibliography   |             |          |             |            | 17  |
| A Fundamentals | of Tensor   | Networks | and Optimal | Truncation | 18  |
A.1 Graphical notation and basic operations . . . . . . . . . . . . . . . . . . . . 18
A.2 The singular value decomposition . . . . . . . . . . . . . . . . . . . . . . . . 19
A.3 Matrix product states and operators . . . . . . . . . . . . . . . . . . . . . . 19
A.4 Gauge freedom and canonical forms . . . . . . . . . . . . . . . . . . . . . . . 20
A.5 Schmidt decomposition and the reduced density matrix . . . . . . . . . . . . 20
A.6 Optimal truncation . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . 21
A.7 The reduced transition matrix and its truncation . . . . . . . . . . . . . . . 21
| B Fisher Zeros | and Dynamical |     | Quantum Phase | Transitions | 22  |
| -------------- | ------------- | --- | ------------- | ----------- | --- |
B.1 Zeros of the partition function . . . . . . . . . . . . . . . . . . . . . . . . . . 22
B.2 The boundary partition function and the Loschmidt amplitude . . . . . . . 22
B.3 Dynamical free energy and dynamical quantum phase transitions . . . . . . 23
B.4 Transfer-matrix realization . . . . . . . . . . . . . . . . . . . . . . . . . . . . 23
2

1 Introduction and Motivation
The study of out-of-equilibrium dynamics in many-body quantum systems remains one of
the most active areas of research in quantum physics, holding the key to understanding
and characterizing fundamental phenomena arising from these complex systems.
A fundamental question in contemporary physics is how a closed quantum system re-
laxes after a sudden perturbation, such as a quantum quench. In generic, non-integrable
systems this relaxation is expected to be governed by the Eigenstate Thermalization Hy-
pothesis (ETH): the system effectively acts as its own heat bath, scrambling information
so that local observables eventually reach thermal equilibrium values that depend only on
macroscopic conserved quantities such as the energy. Integrable systems—which possess
an extensive set of local conserved charges—do not follow this mechanism because they
retain memory of their initial conditions and relax instead to a Generalized Gibbs Ensem-
ble (GGE). Because exact analytical solutions are largely restricted to these integrable
limits, exploring the true mechanisms of thermalization, information scrambling and the
emergence of chaos in non-integrable models relies almost entirely on our ability to track
their long-time dynamics numerically.
Exact numerical simulations are severely limited by the exponential scaling of the
Hilbert space, which makes representing a generic many-body wavefunction intractable.
To circumvent this bottleneck, the community has turned to Tensor Networks (TN), a
framework that efficiently compresses quantum states by using entanglement as the orga-
nizingprinciple. Forone-dimensionalsystemsinequilibrium,MatrixProductStates(MPS)
combined with algorithms such as the Density Matrix Renormalization Group (DMRG)
have become the state-of-the-art approach. Their success comes from the area law of
entanglement: the entanglement entropy across a bipartition scales with the size of the
boundary between the two halves. In 1D, the boundary is just a single point, so the
entropy saturates to a constant instead of scaling with the volume of the chain (which cor-
responds to volume law behaviour). Since the virtual bond dimension required by an MPS
is controlled by precisely this entanglement, a low-entanglement state can be represented
faithfully and efficiently. By restricting the variational manifold to these low-entanglement
states, the MPS ansatz isolates the exact corner of the Hilbert space where equilibrium
physics actually lives.
Extracting dynamics, however, remains a challenge because of the so-called entangle-
ment barrier. During time evolution—for instance, after a quantum quench—the entangle-
mententropyofthestatetypicallygrowslinearlyintime, exhibitingvolume-lawbehaviour.
While highly optimized algorithms such as Time-Evolving Block Decimation (TEBD) or
the Time-Dependent Variational Principle (TDVP) are very successful at short times, this
linear growth forces an exponentially increasing bond dimension to keep the state repre-
sentation faithful, which bottlenecks any standard computation at late times.
To overcome this barrier, recent work has explored a radically different perspective:
transversecontraction. InsteadofevolvingthestateforwardintimeintheusualSchrödinger
picture, the full dynamical evolution of a 1D system is encoded as a 2D tensor network
in which time plays the role of a second spatial dimension. By contracting this network
transversally—alongtheoriginalspatialdirection—onecan,infavourablecases,bypassthe
entanglement barrier and reach late-time observables. In this picture a new quantity ap-
pears naturally: the temporal entanglement entropy, obtained from the standard Schmidt
decomposition taken across a cut in the time domain rather than across the spatial lat-
tice. Its physical interpretation is still an open question, but it is computationally well
defined and it controls the cost of the transverse contraction. nota: debería mencionar la
3

Heisenberg picture, el operator entanglement y su relación con el temporal entanglement?
It is worth mentioning that this cost is dictated by criticality, not by integrability.
When a system is quenched to a conformally invariant critical point, universal data—
such as the central charge of the underlying Conformal Field Theory (CFT)—can be read
directlyfromtheLoschmidtecho, thereturnamplitudeoftheinitialstateafterthequench.
In the transverse picture, this echo is governed by a spatial transfer matrix, and CFT
provides a prediction for its spectrum: at long times, the leading eigenvalues approach
pure phases (their moduli tend to one, with corrections decaying as 1/T, 1/T2,...), so the
transfer matrix becomes effectively unitary. This phenomenon has been called emergent
dual unitarity. Its direct consequence is that the generalized temporal entropies grow only
logarithmically, sothecorrespondingquenchescanbesimulatedwithpolynomialresources.
Carignano and Tagliacozzo [] demonstrated this explicitly for the critical Ising and three-
state Potts models. These are therefore critical-point (CFT) statements: nothing in the
argument requires the model to be exactly solvable.
This raises naturally the central question of the present work: does emergent dual
unitarity, and with it the logarithmic scaling of the generalized temporal entropies, survive
when the model is no longer integrable? To answer it, we study an ANNNI-type model:
a self-dual quantum spin chain that extends the transverse-field Ising model with next-
nearest-neighbour (NNN) interactions controlled by a single coupling p. The model serves
as a deliberate stress-test for the framework. It is non-integrable for any p > 0, yet (as
we verify in Section 2) it remains critical and stays in the Ising universality class over
a wide range of p, so the CFT predictions still have a well-defined target. At the same
time, its longer-range couplings are precisely what makes it expensive for conventional
dynamicalmethods. Inparticular,wepreparethesysteminthefullypolarizedparamagnet
|Ψ ⟩ = |X+⟩⊗N (the ground state of the model in the λ → ∞ limit), quench it to the
0
critical point λ = 1, and use transverse contraction to compute the Loschmidt echo and
the generalized temporal entropies. Comparing them against the Ising CFT prediction lets
us test whether universal signatures and emergent dual-unitary behaviour withstand the
breaking of integrability and the introduction of NNN interactions.
[Thesis outline: Chapter 2 introduces the ANNNI-type model and its equilibrium crit-
icality; Chapter 3 develops the transverse-contraction framework; Chapter 4 details the
exponential-MPO construction; Chapter 5 discusses the non-Hermitian transfer matrix and
dynamical quantum phase transitions; Chapter 6 presents the results.]
2 The ANNNI-type Model
We study a self-dual quantum spin chain that generalizes the transverse-field Ising model
by adding next-nearest-neighbour (NNN) interactions [Alc16]. For the spin-1/2 case its
Hamiltonian reads
H = − X(cid:0) σzσz +pσzσz +λσx+pλσxσx (cid:1) , (1)
i i+1 i i+2 i i i+1
i
where σx,z are the standard Pauli matrices, λ is the transverse field, and p controls the
strength of the longer-range (NNN) coupling.
Setting p = 0 removes the NNN terms and recovers the standard, integrable transverse-
field Ising model, with λ playing the role of the transverse magnetic field. Switching on
p > 0 introduces a second interaction length scale and breaks integrability: the chain
acquires genuine many-body chaos and is no longer exactly solvable.
Breaking integrability raises the following question: does the model remain critical,
and if so, does it stay in the same universality class as the Ising chain? We answer this
4

in equilibrium by examining the ground-state entanglement and comparing it with the
| predictions | of Conformal | Field Theory. |     |     |     |
| ----------- | ------------ | ------------- | --- | --- | --- |
For a critical 1D system, CFT predicts that the ground-state von Neumann entropy of
a block of length l in an open chain of length L scales logarithmically as
|     |     |     | (cid:18)L | lπ(cid:19) |     |
| --- | --- | --- | --------- | ---------- | --- |
c
|     |     | S (l,L) | = ln | sin +k, | (2) |
| --- | --- | ------- | ---- | ------- | --- |
vN
|     |     |     | 6   | π L |     |
| --- | --- | --- | --- | --- | --- |
where c is the central charge of the underlying CFT and k is a non-universal constant.
Evaluating Eq. (2) at the mid-chain cut l = L/2, where sin(π/2) = 1, simplifies it to
|     |     |     | c (cid:18)L(cid:19) | c   |     |
| --- | --- | --- | ------------------- | --- | --- |
lnL+k′.
|     |     | S (L/2) | = ln | +k = | (3) |
| --- | --- | ------- | ---- | ---- | --- |
|     |     | vN      | 6 π  | 6    |     |
Thus, if we run DMRG for a range of system sizes, extract the mid-bond von Neumann
entropy, and plot it against lnL, criticality manifests as a straight line whose slope is
exactly c/6.
Figure 1: Finite-size scaling of the mid-chain von Neumann entropy for the ANNNI-type model at
p=0.5. The linear dependence on lnL is consistent with the CFT prediction of Eq. (3), and from the
| slope we | obtain the central | charge c. |     |     |     |
| -------- | ------------------ | --------- | --- | --- | --- |
As Figure 1 shows, the mid-chain entropy at p = 0.5 follows precisely this linear law,
confirmingthat,evenwithintegrabilitybroken,thesystemstillexhibitsthescaleinvariance
| characteristic | of a CFT. |     |     |     |     |
| -------------- | --------- | --- | --- | --- | --- |
We now ask how the central charge evolves as the NNN coupling is increased. In
his original proposal [Alc16], Alcaraz argued, via finite-size scaling, that for sufficiently
small coupling (p ≲ 1.5) the model stays critical and in the same universality class as the
integrable Ising chain, so the central charge should remain very close to c = 1/2.
To verify this numerically, the central charge was fitted across p ∈ [0,2] on a chain of
N = 300 spins (Figure 2). At the integrable point p = 0 the extracted value is c ≈ 0.502.
The small deviation from the exact c = 1/2 is a well-understood finite-size artifact: on
a finite lattice, the subleading corrections to Eq. (2) have not fully vanished, producing
the systematic offset observed from the thermodynamic limit. As p grows, the fitted
5

Figure 2: Dependence of the fitted central charge c on the frustration parameter p for the Alcaraz model
and a finite chain of N =300 spins. The values remain close to the Ising prediction c=1/2 over the
studied range.
charge drifts gently upward, but the deviation remains small—even at p = 2, it barely
exceeds 0.510. This persistent proximity to c = 1/2 provides strong numerical support for
Alcaraz’s claim [Alc16]: the NNN coupling does not drive the system into a new phase.
The model remains critical and within the Ising universality class, so the CFT framework
is a legitimate tool for describing its dynamics.
Once the equilibrium criticality is established, we move to the dynamical question that
motivates the rest of this work. We quench the system from the polarized initial state
|Ψ ⟩ = |X+⟩⊗N to the critical point λ = 1 and ask whether the universal CFT signatures
0
that organize the equilibrium entanglement also organize the temporal entanglement of
the post-quench dynamics. To address this, we turn to transverse contraction, which gives
us efficient access to the generalized temporal entropies whose CFT predictions we wish to
test.
3 Time Evolution via Transverse Contraction
We implement the transverse contractions with ITransverse.jl, the package developed
by my amazing supervisor Carignano []. The reader interested in the software engineering
is referred there. Here we summarize what the framework computes and how it works.
Transverse contraction translates the time evolution of a D-dimensional system as the
contraction of a (D+1)-dimensional tensor network, the extra dimension being time. For
our 1D chain, this is a 2D network (a Projected Entangled-Pair States, PEPS, up to
Trotter error) of finite bond dimension that encodes the entire dynamics of our system.
The problem then becomes one of contracting it, which is a generally hard task, since
exact PEPS contraction is #P-hard.
In order to visualize this contraction scheme, let us consider a 1D system described by
a Hamiltonian H, and the Loschmidt echo, the return amplitude of an initial state |Ψ ⟩
0
6

after evolving for a time T:
L = ⟨Ψ |e−iHT |Ψ ⟩. (4)
0 0
The ingredients are the initial and final states, written as MPS, and the time-evolution
operator, written as a product of MPO layers. For simplicity, we work with a translation-
invariant, time-independent Hamiltonian that guarantees that all bulk MPO tensors W
are identical (Section 4 explains how to build such an MPO for an NNN model).
In this setup, rather than stacking Trotterized MPO layers on top of the initial MPS—
the usual procedure, which steadily increases the spatial entanglement—we perform a
transverse contraction: we rotate the network by 90◦, so that what were virtual spatial
links become physical temporal links and what were physical spatial links become the vir-
tual temporal bonds. In this rotated picture, the boundary columns are temporal Matrix
Product States (tMPS), ⟨L| and |R⟩, and a bulk column is a temporal Matrix Product
Operator (tMPO), E, which acts as a spatial transfer matrix—which, by translation invari-
ance, is the same fixed operator at every spatial site for different time steps. Contracting
the network along the spatial direction—which is the same as performing time evolution—
then amounts to repeatedly applying E to the boundary tMPS. For an infinite (or long)
chain, this reduces to finding the dominant left and right eigenvectors of E and taking their
overlap, which can be done with the tensor-network Power Method (Section 3.3).
It is worth contrasting the two directions of this network explicitly, because they have
very different spectral character. The rows are the Trotter layers U(δt) = e−iHδt, which
are unitary, and therefore their spectrum lies entirely on the unit circle (the eigenvalues
are pure phases). The column E, on the other hand, is not unitary in general; it is non-
Hermitian, and its eigenvalues are generic complex numbers. The interplay between these
two facts is exactly where the CFT physics enters, as we discuss below and develop in
Section 5.
3.1 Temporal entanglement and the CFT connection
Each application of E increases the bond dimension of the tMPS, so we must necessarily
truncate after every step. Whether this truncation is computationally cheap or not is
determined by the temporal entanglement: the entanglement entropy across a cut in the
time domain of the tMPS. Although its physical meaning is still under investigation, it is
a well-defined computational quantity, and one can study how it scales with the number
of temporal sites N (equivalently, the number of time steps). For many systems, it is
t
observed to obey a volume law, making the contraction no easier than the Schrödinger
approach.
The reason transverse contraction can nevertheless succeed at a critical point is that
the relevant entropies are precisely the quantities that CFT predicts. The key object here
is the Reduced Transition Matrix (RTM), built jointly from the two temporal boundary
states,
T ∝ |R⟩⟨L|, (5)
t
which captures the overlap ⟨L|R⟩ that we actually want to compute. The spectrum of
T defines a family of generalized temporal entropies, and the central result we rely on
t
is that, for a quench to a conformally invariant critical point, these generalized entropies
are exactly the entropies computed in CFT, which grow only logarithmically with time [].
This means the bond dimension needed to represent the tMPS grows only polynomially,
and a truncation that targets the RTM spectrum keeps the entire transverse contraction
efficient. In other words, the efficiency of the RTM approach is not primarily algorithmic,
7

but physical: it is a consequence of its spectrum being fixed by the CFT at the critical
point.
| 3.2 Truncation |     | based | on  | RDM | and RTM |     |     |     |     |     |
| -------------- | --- | ----- | --- | --- | ------- | --- | --- | --- | --- | --- |
Nota: explicar bien a fondo las canonical forms y métodos de truncación básica en RDM
en el APÉNDICE
We have two complementary truncation schemes at our disposal, based respectively on
the Reduced Density Matrix (RDM) and the Reduced Transition Matrix (RTM). The first
is the conventional truncation available in most tensor-network libraries; the second is the
gauge-invariant scheme that ITransverse.jl implements for the transverse picture. They
differ in two respects that turn out to be linked—gauge invariance and bond-dimension
| cost—and | each | has | a regime | in  | which it | is preferable. |     |     |     |     |
| -------- | ---- | --- | -------- | --- | -------- | -------------- | --- | --- | --- | --- |
The RDM truncation acts on a single boundary state: one brings the MPS to mixed
canonical form and discards the smallest singular values of the orthogonality center, which
istheoptimallow-rankapproximationofthatstateinthe2-norm. Inthetransversepicture,
however,thetransfermatricesarenon-Hermitian,andtheRDMtruncationbecomesgauge-
dependent. In order to visualize why, consider the overlap ⟨L|R|L|R⟩ (or ⟨L|E|R⟩ for an
operator) which encodes the physically relevant content of the network. Thanks to the
gauge freedom within our tensor network, we can insert an identity XX−1 = I between
adjacent columns while keeping the exact global contraction invariant,
|     |     |     |           |     |            |     |     | D (cid:12)         | E   |     |
| --- | --- | --- | --------- | --- | ---------- | --- | --- | ------------------ | --- | --- |
|     |     |     |           |     | ⟨L|XX−1|R⟩ |     |     | L˜|R˜(cid:12)L˜|R˜ |     |     |
|     |     |     | ⟨L|R|L|R⟩ |     | =          |     | =   |                    | .   | (6) |
(cid:12)
The boundary RDM, by contrast, is not invariant: for the gauged left state,
|     |     |     |     |     | (cid:12)        | ED (cid:12)  |     |     |     |     |
| --- | --- | --- | --- | --- | --------------- | ------------ | --- | --- | --- | --- |
|     |     |     |     |     | ρ˜ = (cid:12)L˜ | L˜(cid:12) = | X†ρ | X,  |     | (7) |
|     |     |     |     |     | L (cid:12)      | (cid:12)     |     | L   |     |     |
isingeneralnotunitary(X†X
| andsinceX |     |     |     |     |     | ̸= I),itsspectrumcanbesignificantlydistorted. |     |     |     |     |
| --------- | --- | --- | --- | --- | --- | --------------------------------------------- | --- | --- | --- | --- |
The isolated RDM spectrum is therefore not a canonical, gauge-invariant quantity: if we
perform a truncation by discarding the smallest singular values of these isolated RDMs, it
is possible that the singular values that carry the critical information required to preserve
theglobalscalarproduct⟨L|R⟩mightbepushedintothetailsofthespectrumbythegauge
transformationandsubsequentlythrownaway, fundamentallybreakingthecalculation. As
long as we are in a fixed canonical form, the RDM truncation remains well defined and
usable. While this is the standard scheme, it does not target the overlap directly, which is
| the conceptual |     | reason | to look | for | an alternative. |     |     |     |     |     |
| -------------- | --- | ------ | ------- | --- | --------------- | --- | --- | --- | --- | --- |
The RTM approach removes the gauge ambiguity entirely. Under the same gauge,
|     |     |     |     | (cid:12) ED | (cid:12)   |            |     |        |     |     |
| --- | --- | --- | --- | ----------- | ---------- | ---------- | --- | ------ | --- | --- |
|     |     |     | T˜  | (cid:12)R˜  | L˜(cid:12) | X−1|R⟩⟨L|X |     | X−1TX, |     |     |
|     |     |     |     | =           | =          |            |     | =      |     | (8) |
|     |     |     |     | (cid:12)    | (cid:12)   |            |     |        |     |     |
T˜
which is a similarity transformation: and T share the same eigenvalues, so the RTM
spectrum is gauge-invariant. The RTM is therefore both the object CFT predicts and the
| only gauge-stable |     | thing | to  | truncate | on. |     |     |     |     |     |
| ----------------- | --- | ----- | --- | -------- | --- | --- | --- | --- | --- | --- |
In practice, neither scheme requires computing the eigenvalues of the exponentially
large global matrix: thanks to canonical forms, it is possible to reduce the problem to
diagonalising a local χ×χ environment matrix, M . For the Hermitian RDM, a bipartite
A
cut gives ρ = U M U†, with U being unitary (built from left-orthogonal tensors), so
|     | B   | B   | A B |     | B   |     |     |     |     |     |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
eig(ρ ) = eig(M ), and one can truncate rigorously and straightforwardly on the eigen-
| B   |     | A   |     |     |     |     |     |     |     |     |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
values of the smaller matrix. For the RTM, built from two different states, the same cut
8

gives τ = U M V , with U ̸= V†. This is no longer a similarity transformation as
B B A B B B
before, so eig(τ ) ̸= eig(M ). Fortunately, U and V are isometries, so they preserve
B A B B
singular values: sv(τ ) = sv(M ), allowing us to truncate on the singular values of the
B A
small environment instead.
This last step requires special care, though. The singular values of a matrix M are
q
always related to its eigenvalues via sv(M) = eig(MM†). For a Hermitian positive-
definiteRDM,sincethesingularvaluescoincidewiththeeigenvalues,droppingthesmallest
singular values is mathematically equivalent to dropping the smallest eigenvalues. For the
non-Hermitian RTM, however, the eigenvalues can be complex or negative, and the trace
TrT = P λ —which gives us the physical overlap we are interested in computing—is no
t i i
longer in one-to-one correspondence with the singular-value magnitudes. To see why, let
us consider the toy example of an RTM with eigenvalues +1 and −1: while their singular
values (their moduli) are both large (1) and would therefore be preserved by the Singular
ValueDecomposition(SVD)algorithm, theyexactlycanceloutwhenevaluatingtheactual
physical quantity of interest. RTM truncation therefore lacks the rigorous mathematical
boundsofRDMtruncationandwouldlikelyfailforagenericrandomnon-Hermitianmatrix.
The reason why it nevertheless works in practice is physical: because the transfer matrices
related to time evolution have highly constrained, physical structures, these cancellations
are in practice harmless.
OneadvantageoftheRTMapproachisthatitisconsiderablycheapertocomputethan
the RDM one. We can understand better why from what each scheme tries to preserve.
The eigenvalues of the boundary RDM ρ = Tr |R⟩⟨R| are the squared Schmidt coef-
R A¯
ficients of |R⟩ across the cut, so its truncation keeps the directions in which |R⟩ carries
the most weight, regardless of whether those directions overlap with ⟨L|. It thus faithfully
represents the entire state |R⟩, including the component orthogonal to ⟨L| that will con-
tribute nothing once the overlap ⟨L|R⟩ is taken. The RTM, T ∝ |R⟩⟨L|, instead weights
each direction by its contribution to the overlap itself: a direction along which |R⟩ is large
but which is orthogonal to the corresponding part of ⟨L| carries a small singular value and
is thus discarded. By retaining only the part of each boundary that actually contributes
to ⟨L|R|L|R⟩, the RTM needs fewer Schmidt values—which translates into a lower bond
dimension—than the RDM, which pays to store information the final overlap throws away.
Formally, the RDM compression is bounded by the standard state entanglement, whereas
the RTM compression is bounded by the smaller operator-space entanglement [].
In summary, the two schemes are complementary and we use both. The RTM trunca-
tion is gauge-invariant and substantially cheaper in bond dimension, and it is our default
for production runs; but it rests on a singular-value heuristic with no strict bound on the
overlap, so its convergence can be less smooth and occasionally unstable, and it must be
monitoredwithcare. TheRDMtruncationisgauge-dependentandmoreexpensive, butas
the truncation of a Hermitian, positive operator, it inherits a rigorous variational bound,
so it converges smoothly and monotonically and yields more robust results. We therefore
fall back to the RDM scheme whenever the RTM output is noisy or difficult to interpret,
accepting the higher bond dimension in exchange for a cleaner and more reliable answer.
3.3 The Power Method
Ultimately, we are interested in computing the full transverse contraction of the (infinite)
2Dtensornetworkencodingthedynamicsofoursystem. Forspatiallytranslation-invariant
systems, this task reduces to finding the dominant left and right eigenvectors of a single
spatial transfer matrix E and computing their overlap ⟨L|R|L|R⟩. To achieve this, we will
9

use the tensor-network version of the Power Method.
Starting from initial guesses for the temporal boundary states (tMPS) ⟨L | and |R ⟩,
0 0
we repeatedly absorb at each step i one column of the network—corresponding to our
transfer matrix—into each boundary,
D (cid:12) (cid:12) E
L˜ (cid:12) = ⟨L | E, (cid:12)R˜ = E |R ⟩. (9)
i+1(cid:12) i (cid:12) i+1 i
Each application multiplies the bond dimension, so we truncate (on the RTM or RDM)
after every step. In addition, because E is non-unitary, we also renormalize at every step
to avoid numerical instabilities. The iteration process stops once the generalized temporal
entropy stabilizes between consecutive steps, ∆S < ε, with ε taken small. Because E
is non-Hermitian, its left and right eigenvectors differ, so the two boundaries must be
evolved independently—which is achieved by the left–right power method implementation
in ITransverse.jl, appropriate to the asymmetric MPO of our model (Section 4).
It is important to consider some caveats regarding the choice of the truncation scheme
used for the power method. Ideally, we would truncate using the full remaining network,
e.g. the global transition matrix ⟨L|E···E|R⟩. But the power method grows the net-
work iteratively column by column, so whether we are extracting the dominant singular
values from the local RDM of |R˜ ⟩ or from the joint local RTM governing the overlap
i+1
⟨L˜ |R˜ ⟩, the algorithm optimizes at each step the boundary states based exclusively
i+1 i+1
on the current local environment alone. It completely ignores the remaining E columns
that have yet to be applied. One might worry that information which looks irrelevant
now (small local singular values) but is relevant later for the full contraction is discarded
prematurely. But for homogeneous problems—such as the Loschmidt echo of a translation-
invariant chain—this is not a problem: every bulk column is the same E, so if relevant
weight is dropped at step i, the next application of the identical E regenerates it, and
the method converges ultimately to the true dominant eigenspace (provided the system is
gapped; at a closing gap it converges to the degenerate subspace, a point we return to in
Section 5).
3.4 Generalized temporal entropies
There is an important relationship between the entanglement entropy of a state and the
bond dimension required to faithfully represent it with a tensor network. In standard
Schrödinger time evolution, the complexity is governed by the spatial entanglement of the
state, which typically grows linearly with time (volume law). In the transverse picture, on
theotherhand,theobjectsofinterestarenowtemporalboundarystates,whosecomplexity
is encoded in the spectrum of the RTM T ∝ |R⟩⟨L|. This leads us to the concept of
t
generalized temporal entropies, first introduced in high-energy physics and holography to
describe transition matrices between different quantum states—which correspond to our
left, ⟨L|, and right, |R⟩, boundary states.
From the RTM T , we can define the generalized von Neumann temporal entropy as:
t
Sgen(t) = −Tr(T logT ). (10)
1 t t
This definition presents several practical challenges, primarily inherited from the math-
ematical nature of the RTM. Because T is non-Hermitian, its eigenvalues—and hence
t
Sgen—can be complex, requiring its physical interpretation to be handled with care. Fur-
1
thermore, evaluating the matrix logarithm requires the exact eigenvalues of T , which is
t
an exponentially large matrix. As we discussed in the RTM and RDM truncations sec-
tion, our local tensor network algorithm only gives us access to the singular values of the
10

smallerenvironment. Butbecausethereisnoexactone-to-onecorrespondencebetweenthe
singular values and the eigenvalues of a non-Hermitian operator—unlike the well-behaved
Hermitian RDM—extracting the exact von Neumann entropy becomes a very complicated
task (although it is possible to do it, we refer to ITransverse paper).
We therefore work instead with higher-order generalized entropies, such as the general-
ized Rényi-2 entropy,
Sgen(t) = −logTr(T2). (11)
2 t
which needs neither a diagonalization nor a logarithm: Tr(T2) is obtained by taking two
t
copies of the RTM and contracting them, making the calculation efficient and numerically
robust. These generalized entropies are of central interest because they dictate the cost
of the transverse contraction: whereas spatial entanglement always grows linearly in time,
the generalized temporal entropy can behave far more favourably—and, as discussed in
Section 3, CFT predicts logarithmic growth at a critical point.
Operationally, once the power method has converged to ⟨L| and |R⟩ for evolution
up to a final time T, we evaluate Sgen across every bipartition of the temporal chain:
2
each internal bond corresponds to a temporal cut t ∈ (0,T). A dedicated routine in
ITransverse.jl builds the two RTM copies and contracts them at each link. Because the
transition matrices are non-Hermitian, the resulting entropy is generally complex, so we
plot the real and imaginary parts separately against the (normalized) temporal cut. This
profiletellsushowthecomplexityofthecontractionevolveswithtimeand—bycomparison
withtheCFTprediction(Section??)—whethertheuniversalIsingbehaviourpersistsonce
integrability is broken.
4 MPO Representation of e−iHt
We now construct an MPO for the time-evolution operator U(δt) = e−iHδt of the ANNNI-
type model that is compatible with transverse contraction. As already emphasized in
Section 3, the MPO should have translation invariance: this way, when the network is
rotated and the columns of the time-evolution MPO become the spatial transfer matrices,
the whole 2D network is simply the infinite repetition of a single identical column. Only
thendoesthetransversecontractionreducetoaneigenproblemforjustonetransfermatrix,
giving easy access to the thermodynamic limit. Any site dependence in the MPO would
require a different transfer matrix at every site and destroy the efficiency and simplicity of
the method—in this case, of our version of the power method.
For standard nearest-neighbour (NN) models, such as the transverse-field Ising model,
generating this translation-invariant MPO is straightforward. The usual approach is to
apply a Trotter–Suzuki decomposition in small time steps, δt, so that the full evolution
is written as a product of short-time operators, U(t) = QU(δt). Each U(δt) is then
split into alternating layers of even and odd two-site gates, forming the familiar “brick-
wall” circuit. This works well because every gate acts on adjacent physical sites: pairs
of decomposed tensors can be contracted analytically and reshaped into a uniform one-
dimensional sequence of identical MPO tensors (see ITransverse paper for details).
The NNN coupling of our model, by contrast, breaks this simple picture. An interac-
tion connecting site i directly to site i+2 does not live on a single physical bond, so it
cannot be placed into the same local two-site pattern as in the NN case. A Trotterized con-
struction would then have to either group three sites into one enlarged tensor or introduce
SWAP gates to bring the non-adjacent spins next to each other. Both options make the
circuitmorecumbersomeand,moreimportantlyforus,destroythecleanuniformstructure
required to build an efficient translation-invariant MPO.
11

| 4.1 The | Hamiltonian | as  | a Finite-State |     | Machine |     |     |     |     |
| ------- | ----------- | --- | -------------- | --- | ------- | --- | --- | --- | --- |
To keep a perfectly uniform, translation-invariant MPO for the NNN case, we must aban-
don spatial Trotterization and use a different strategy: we encode the full Hamiltonian
directly as a block-upper-triangular MPO. In this representation, the virtual indices of the
MPO behave like the internal states of a finite-state machine that records (into a mem-
ory) whether an interaction has been started, propagated, or terminated. For a generic
extensive Hamiltonian, the local tensor can be written in the following block form
|     |     |     |     |     |    |    |     |     |      |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | ---- |
|     |     |     |     |     | I   | C D |     |     |      |
|     |     |     |     |     | 0  | B, |     |     |      |
|     |     |     |     | W H | =   | A   |     |     | (12) |
|     |     |     |     |     |    |    |     |     |      |
|     |     |     |     |     | 0   | 0 I |     |     |      |
where D contains the local on-site terms, C and B initiate and terminate interactions,
respectively, and A propagates an interaction across the intermediate sites.
FortheANNNI-typemodel,theseblocksmustaccommodatethetransversefield(−λσx),
the two NN interactions (−σzσz and −pλσxσx), and the NNN coupling (−pσzσz). The
| resulting | local MPO | tensor | has bond |     | dimension | D = | 5 and | reads |     |
| --------- | --------- | ------ | -------- | --- | --------- | --- | ----- | ----- | --- |
w
|     |     |     |    | −σz | −pλσx | −pσz | −λσx |     |     |
| --- | --- | --- | --- | --- | ----- | ---- | ----- | --- | --- |
I
|     |     |     | 0  | 0   | 0   | 0   |     | σz  |     |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |

|     |     |     |     |     |     |     |     |    |      |
| --- | --- | --- | ---- | --- | --- | --- | --- | --- | ---- |
|     |     | W   | = 0 | 0   | 0   | 0   | σx  | .  | (13) |
|     |     | H   |     |     |     |     |     |    |      |
|     |     |     | 0   |     |     |     |     |    |      |
|     |     |     |      | I   | 0   | 0   |     | 0   |      |
|     |     |     |     |     |     |     |     |    |      |
|     |     |     | 0    | 0   | 0   | 0   |     | I   |      |
Here D = −λσx contains the on-site transverse field. The first row contains the inter-
action initiators, C = (−σz, −pλσx, −pσz); the last column contains the terminators,
| (σz, | σx, 0)T; |     |     |     |     |     |     |     |     |
| ---- | -------- | --- | --- | --- | --- | --- | --- | --- | --- |
B = and the internal block A stores the memory channels that propagate
| unfinished | (NNN) | interactions | across |     | the virtual | bond. |     |     |     |
| ---------- | ----- | ------------ | ------ | --- | ----------- | ----- | --- | --- | --- |
It is useful to read this MPO as a directional finite-state machine (FSM) that scans
the chain from left to right, using the virtual bond as its internal memory. The machine
starts in the identity state (the first row). From there, it can either apply an on-site
operator through D and jump directly to the final completed state, or it can open a multi-
site interaction through the initiator block C and enter the memory block A. For a NN
interaction, the machine simply waits one site and terminates through B. The NNN term
is where the FSM construction becomes especially useful: initiating it puts the machine
intoadedicatedmemorystate(thefourthrowofEq.(13)),whichappliesanidentityatthe
σz
intermediate site i+1 and carries the pending operator untouched until the interaction
terminates correctly at site i+2. Once an interaction is finished, the machine moves into
the final absorbing state, where it applies identities to the rest of the chain.
Because the matrix is strictly upper-triangular, contracting the MPO from left to right
can only move forward through this sequence; it can never run an interaction backwards.
As we have seen, this perfectly accounts for the long-range NNN terms while preserving
translation invariance along the MPO tensors. On the other hand, the required bond
dimension is simply the number of internal machine states: the initial state, the final
absorbing state, and one memory channel for each propagating interaction, giving D = 5
W
| for the | ANNNI Hamiltonian. |     |     |     |     |     |     |     |     |
| ------- | ------------------ | --- | --- | --- | --- | --- | --- | --- | --- |
A consequence of this directional structure, important for what follows, is that the
MPO explicitly breaks left–right parity. Because the finite-state machine reads the chain
only in one direction, the block matrix is inherently asymmetric in its virtual bond space.
When the MPO is rotated by 90◦ to form the temporal MPO (tMPO), the corresponding
12

spatialtransfermatrixisthereforenon-Hermitianandstructurallyasymmetric. Symmetric
MPOs admit useful shortcuts—for instance, the right temporal state can be obtained by
transposing the left one—but generic asymmetric MPOs like ours require independent left
and right boundary states. This is precisely the setting handled by the left–right power
method implemented in ITransverse.jl, and we return to the physics of the resulting
non-Hermitian transfer matrix in Section 5.
4.2 Exponentiating the Finite-State Machine
check this section again
Once the Hamiltonian is encoded in the finite-state-machine MPO W , the remaining
H
task is to turn it into a time-evolution MPO W(τ) ≈ eτH, with τ = −iδt. This step, which
involves exponentiating W , is not straightforward, though. Because the block is upper
H
triangular,itspowersmixthedifferentblockstogether. Forinstance,ifwewanttocompute
(τW )2,crosstermsstartappearingsuchasCD,BC,AD,ogetherwiththeircounterparts
H
DC, CB, DA. Each of these must appear with exactly the coefficient prescribed by the
Taylor series of eτWH, or the resulting operator stops matching that series at the order we
are targeting, and the truncation error no longer scales as the expected power of δt.
There is also a deeper numerical obstacle. At first, one might think of truncating the
Taylor series eτH = I+τH + τ2 H2 +··· at some order; but this cannot be carried out
2!
directly in the thermodynamic limit. Each power Hn is a sum of O(Nn) local terms,
so successive orders scale with different powers of the chain length N, and the resulting
state cannot be normalized as N → ∞ [VDHMV24]. Any acceptable construction must
therefore be size-extensive: applied to a normalizable state, it has to return a normalizable
state regardless of how long the chain is.
The key idea to get around both problems is to organize the expansion not by powers
of H, but by how its local terms overlap along the chain. Terms acting on disjoint regions
commute, and their product factorizes trivially. These disconnected contributions are
harmless and can be kept to all orders at no cost. The true difficulty lies in the connected
contributions, where two or more local terms might act on overlapping sites and need
not commute. For instance, an on-site field term D may act on the same site where
an interaction is being initiated by C. Both orderings, CD and DC, are then required:
H2 = (P h )2 sums over every ordered pair of local terms regardless of whether they
i i
commute, and dropping either ordering would in any case break the Hermiticity that H2
must inherit from H, since individually (CD)† = DC ≠ CD. A suitable time-evolution
MPO must therefore keep all disconnected terms exactly, capture the connected clusters
correctly up to some controlled size, and push the truncation error to a chosen, known
order in δt.
The simplest such operator keeps every disconnected term and discards all overlaps.
It has an exact and compact MPO—the WI operator of Zaletel et al. []—built directly
from the blocks A,B,C,D of W by distributing the time step across the initiators C, the
H
terminators B, and the on-site block D (we make this split precise, and explain why it is
builtthewayitis,oncewereachthenumericalimplementationbelow). Itischeap,keeping
the virtual bond at 1+χ, where χ is the number of memory channels of the Hamiltonian
(χ = 3 for our model). However, since it drops every overlapping term, it is only first order.
Over many time steps, the non-unitarity of this approximation can accumulate and make
conserved quantities and observables slowly drift away, producing the so-called “Trotter
leak”.
Wecandoconsiderablybetteratthesamebonddimension. Theimprovement,alsodue
to Zaletel et al., is to additionally keep the terms that overlap on a single site, capturing in
13

particular the on-site contributions to all orders. The obstacle is purely algebraic: writing
out the correct symmetric combinations of overlapping blocks by hand quickly becomes a
cumbersome task. Their clever solution is to encode the four blocks in the dynamics of
two auxiliary hard-core bosons: one exponentiates a small effective Hamiltonian built from
A,B,C,D acting on these bosons, and then reads off the new MPO blocks by projecting
the result onto the boson occupations. Because a hard-core boson can be occupied at most
once, theexponentialtruncatesonitsownandproducesexactlythesymmetriccrossterms
required. The resulting operator WII keeps the on-site term exactly (eτD) and has a bond
dimension no larger than WI, which makes it the natural improved first-order choice.
It is important to note a subtlety regarding the dynamics of our model. Because
WII captures all single-site overlaps exactly, it behaves as a genuine second-order scheme
whenever the Hamiltonian is strictly nearest-neighbour: in that case, two local interaction
terms cannot overlap on more than one site. Our ANNNI-type model, however, is not
nearest-neighbour. The NNN bridge spans three sites, so two such terms can overlap
over a larger region. In particular, a single site can simultaneously host a transverse
field, the beginning of a new interaction, and the propagation of a previous one. These
longer connected clusters are precisely the terms that WII does not reproduce. For this
reason, this construction applied to our Hamiltonian is only first order, which would force
impractically small time steps to keep the accumulated error under control over a long
evolution.
To recover true second-order accuracy independently of the interaction range, we use
the systematic higher-order construction of Van Damme et al. [VDHMV24], which gener-
alizes this auxiliary-space idea to any order following this exact connected/disconnected
logic. Their order-n MPO reproduces all connected clusters of overlapping terms up to a
controlled size, together with all disconnected products of such clusters, and is accurate
to order n in the time step, with a one-step error O(δtn+1). Its first-order MPO coincides
exactlywithWII andhandlesoverlapsconfinedtoasinglesite. Theoneweneedandusein
our calculations is the second-order MPO, denoted “VD2”, which additionally reproduces
the longer overlaps generated by the next-nearest-neighbour terms.
This accuracy comes at a modest price in bond dimension. While WI and WII keep
the virtual bond at 1+χ, the VD2 operator enlarges it to 1+χ+χ2 (here 1+3+9 = 13):
the extra χ2 block stores the doubled memory channel needed to follow two interactions
that overlap at once. We do not reproduce the explicit form of VD2—the closed-form
expressions for its blocks are tabulated in Appendix A of Ref. [VDHMV24], and can also
befoundintheAppendixbelow—andsimplyuseitasablackboxthat,foragivenmemory
dimension χ, returns the time-evolution MPO of the prescribed order. We can now also
make precise the time-step split across C and B that we mentioned earlier: the same
convention underlies all three schemes alike, and it is what keeps them numerically stable.
√
The complex time step τ = −iδt is not divided evenly. We assign a real factor t = δt
C
to the initiator block C and put the complex phase in t = τ/t on the terminator block
B C
B, so that the product C ⊗B still carries the correct factor −iδt while the intermediate
memory channels remain numerically well behaved. The structural fact we carry forward
is the enlarged VD2 bond: under the 90◦ rotation of Section 3, the spatial virtual bond of
U(δt) becomes the physical dimension of the temporal sites, so the tMPS inherits a larger
physical dimension for VD2 than for WII. In exchange, being genuinely second order,
VD2 lets us reach long times with comparatively large steps (δt = 0.05) while keeping the
Trottererrorcontrolled, allowingustoaccesslongertemporalnetworkswiththetransverse
contraction framework.
14

5 The Non-Hermitian Transfer Matrix and Dynamical Phase Transitions
OneofthemostimportantobjectsinthetransversepictureisthespatialtransfermatrixE,
a column of the rotated network. Unlike the Hamiltonian, or the time-evolution operator
U(δt), this matrix is generally non-Hermitian—and, due to the asymmetry of the MPO
(Section4), itisalsonon-normal(EE† ̸= E†E). Non-Hermitianoperatorsbehavedifferently
from the Hermitian operators of standard quantum mechanics, and these differences are
important both for the physics we want to extract and for the numerical methods we use
below.
5.1 Left and right eigenvectors, and biorthogonality
A Hermitian operator has a single orthonormal eigenbasis. A non-Hermitian E has two:
right eigenvectors, defined by E|R ⟩ = µ |R ⟩, and left eigenvectors, by ⟨L |E = µ ⟨L |,
i i i i i i
and in general |L ⟩ ̸= |R ⟩. These sets are not individually orthonormal; instead they
i i
form a biorthogonal system, ⟨L |R |L |R ⟩ ∝ δ . This has two direct consequences for
i j i j ij
the calculation. First, the dominant left and right temporal states must be obtained
independently—therightstatecannotbeobtainedbyconjugatingtheleftoneashappends
inthesymmetriccase. Second,thephysicallymeaningfulcontractionisthebilinearoverlap
⟨L|R|L|R⟩ that encodes the network, not the usual sesquilinear overlap ⟨L∗|R⟩. Therefore,
overlaps must be taken without complex conjugation.
Since the eigenvalues µ are generally complex, we order them by their modulus. The
i
eigenvalue with the largest modulus, µ , controls the long-chain contraction and therefore
0
determines the Loschmidt rate in the thermodynamic limit:
1
ℓ(T) = − lim log|L(T)| = −log|µ (T)|. (14)
0
N→∞ N
For a finite open chain, this relation receives O(1/N) boundary corrections. Thus, when
comparing−log|µ |withafinite-sizeLoschmidtechodividedbyN,theseboundaryeffects
0
have to be taken into account.
5.2 Eigenvalues versus singular values
For Hermitian positive operators, such as the RDMs used in conventional tensor-network
methods, eigenvalues and singular values carry essentially the same information: the eigen-
values are real, positive, and coincide with the squared singular values. This is what makes
the usual truncation procedure variationally controlled. For the non-Hermitian transfer
matrix E, and for the associated RTM, this connection is lost. Its trace, TrE = P µ , is a
i i
sum of complex eigenvalues, so different contributions can cancel each other even when the
singular values are large (as in the {+1,−1} example of Section 3). As a result, truncating
by singular values does not provide the same direct control over the physical contraction.
The situation is further complicated by non-normality: finite applications of a non-normal
matrix can show transient amplification that is not visible from the eigenvalues alone. In
practice these effects remain manageable for the transfer matrices we study, but they ex-
plainwhyRTMtruncation, althoughmoreefficient, islessrigorouslycontrolledthanRDM
truncation.
5.3 Emergent dual unitarity as a spectral statement
Emergent dual unitarity is most naturally seen in the spectrum of the spatial transfer
matrix E. The important point is to separate two pieces of information carried by each
15

eigenvalue: its modulus, which controls growth or decay under repeated contractions, and
| its phase, which | carries | the oscillatory | dynamics. |     |     |
| ---------------- | ------- | --------------- | --------- | --- | --- |
For a unitary operator U, all eigenvalues lie on the unit circle. They can be written
eiθj,
as µ j = so |µ j | = 1 for every j. In this case there is no hierarchy in modulus: no
eigenvalue is larger than the others, and the spectrum is distinguished only by its phases.
This is true for the rows of our two-dimensional network, namely the time-evolution layers
e−iHδt,
| U(δt) = | which | are unitary by | construction. |     |     |
| ------- | ----- | -------------- | ------------- | --- | --- |
The column transfer matrix E is different. At finite time it is not generally unitary,
so its eigenvalues need not stay on the unit circle. They typically move into the complex
plane and develop a hierarchy in modulus. The eigenvalue with the largest modulus, µ ,
0
thendominatesthe contraction, while the subleading eigenvaluesaresuppressedbypowers
of |µ /µ | < 1. This is why the Loschmidt rate is fixed by ℓ(T) = −log|µ (T)| (Eq. (14)),
1 0 0
and why the gap ratio |µ |/|µ | also controls the convergence of the power method.
1 0
After a quench to a critical point, conformal invariance predicts that this hierarchy
gradually disappears as the temporal extent T grows. The leading eigenvalues move back
towarda common circlein the complex plane: their moduli become degenerate, whiletheir
| phases remain | distinct. | In other words, |                 |     |      |
| ------------- | --------- | --------------- | --------------- | --- | ---- |
|               |           | |µ (T)|         | a b             |     |      |
|               |           | i               | i i             |     |      |
|               |           | =               | 1− − −··· −−−−→ | 1,  | (15) |
|               |           | |µ (T)|         | T T2            |     |      |
|               |           | 0               | T→∞             |     |      |
with corrections that decay as integer powers of 1/T. A band of eigenvalues with equal
moduli and different phases is the spectral signature of a unitary operator, up to an overall
scale. This is the sense in which the spatial transfer matrix becomes effectively unitary at
long times. Importantly, this prediction comes from conformal invariance at the critical
| point, not from | integrability. |     |     |     |     |
| --------------- | -------------- | --- | --- | --- | --- |
This also explains why we monitor the leading band of eigenvalues rather than only µ .
0
The approach to dual unitarity shows up in two complementary ways: the leading moduli
collapse toward a common value, equivalently |µ |/|µ | → 1, and the corresponding phases
1 0
continue to spread around the unit circle as T changes. Once several eigenvalues are nearly
degenerate in modulus, a single dominant eigenvalue is no longer clearly separated, so a
one-vector power method becomes unreliable. One must instead track the leading band as
a whole.
Thesamenear-degeneracyisalsothemechanismbehindthedynamicalquantumphase
transitionsdiscussedbelow. AsT changes,twoleadingeigenvaluescanexchangethelargest
modulus, making µ (T), and therefore the rate ℓ(T), non-analytic. Thus the collapse of
0
theleadingmoduliisboththesignatureofemergentdualunitarityandthespectralsetting
| in which DQPTs       | appear. |         |     |     |     |
| -------------------- | ------- | ------- | --- | --- | --- |
| 6 Conclusions        | and     | Outlook |     |     |     |
| yet to be written... |         |         |     |     |     |
16

Bibliography
[Alc16] F. C. Alcaraz. Universal behavior of the shannon mutual information in
nonintegrable self-dual quantum chains. Phys. Rev. B, 94:115116, 2016.
[CRT24] Stefano Carignano, Carlos Ramos Marimón, and Luca Tagliacozzo. On tem-
poralentropyandthecomplexityofcomputingtheexpectationvalueoflocal
| operators after | a quench. | Phys. Rev. | Research, 6:033021, | 2024. |
| --------------- | --------- | ---------- | ------------------- | ----- |
[HPK13] Markus Heyl, Anatoli Polkovnikov, and Stefan Kehrein. Dynamical quan-
tum phase transitions in the transverse-field ising model. Phys. Rev. Lett.,
| 110:135704, | Mar 2013. |     |     |     |
| ----------- | --------- | --- | --- | --- |
[VDHMV24] MaartenVanDamme,JuthoHaegeman,IanMcCulloch,andLaurensVander-
straeten. Efficient higher-order matrix product operators for time evolution.
| SciPost Phys., | 17:135, | 2024. |     |     |
| -------------- | ------- | ----- | --- | --- |
17

A Fundamentals of Tensor Networks and Optimal Truncation
make it more unique, personal, and less generic
The simulation of a quantum many-body system is fundamentally bottlenecked by
the “curse of dimensionality”: for N spin-1/2 particles, the Hilbert space has dimension
2N, meaning a completely general wavefunction requires an exponentially large number
of coefficients. However, nature is remarkably efficient, and Tensor Networks exploit this
fact. Insteadofattemptingtheimpossibletaskofrepresentingarbitrarystates, theytarget
the physically relevant corner of the Hilbert space, where states typically exhibit very low
entanglement. For example, the ground states of local gapped Hamiltonians obey what
is called an “area law”: the entanglement across a cut grows only with the boundary of
the subsystem, not its volume. This allows us to replace one enormous rank-N tensor—
representing the coefficients of the many-body wavefunction—with a network of smaller,
tractabletensorsconnectedbyinternalbondsofcontrolleddimension,givenbyD. Thegoal
of this appendix is to introduce the foundational tools that make this mathematical ansatz
so powerful: graphical notation, singular value decompositions, matrix product states,
canonical forms, optimal truncation, and finally, the non-Hermitian extension required
to overcome the entanglement barrier via the transverse contraction scheme detailed in
Section 3.
A.1 Graphical notation and basic operations
Tensor-network calculations quickly become difficult to read if every tensor is written only
in index notation. Even a moderate contraction, such as
X
C = A B D , (16)
αβγδ αβij jkγ kiδ
i,j,k
contains enough indices that the structure of the operation is no longer immediately clear.
The expression is correct, but it is hard to see at a glance which tensors are connected,
which indices remain open, and which contraction order would be efficient.
For this reason, the standard diagrammatic notation for tensor networks—originally
introduced by Penrose—is incredibly useful. A tensor is then represented by a shape with
one leg for each index. A scalar has no legs, a vector has one, a matrix has two, and a
rank-R tensor T has R legs. In this notation, the rank of an object is visible directly
i1i2...iR
from the picture, rather than being hidden in a long list of indices.
Contractions are represented by joining legs. If two tensors share an index and we sum
overit, thecorrespondinglegsareconnectedbyaninternalline. Thelegsthatremainopen
become the free indices of the resulting tensor. This gives a direct visual representation of
the same algebraic operation as matrix multiplication, but it also makes the structure of
larger tensor networks much easier to follow.
The diagram is not only a convenient shorthand. It also helps identify the contraction
path, this is, the order in which the tensors should be contracted for efficiency. Different
orders give the same final tensor, but their computational costs can be very different: a
poor choice may create a very large intermediate object, while a better one avoids it. For
this reason, choosing an efficient contraction sequence is an essential part of developing
tensor-network algorithms.
A second basic operation is reshaping. Linear-algebra routines, such as eigensolvers
and singular value decompositions, act on matrices, so higher-rank tensors must often be
viewed as matrices by grouping several indices into composite ones. For example, a rank-3
tensor A can be reshaped into a matrix M by treating (i,α) as a single index. We
i,α,β (iα),β
18

simplychangethewaywelabelthesameentries,nottheircontent. Thissimpleoperationis
what allows us to apply matrix tools, especially the singular value decomposition, directly
| inside tensor | networks. |       |               |     |     |     |     |     |     |
| ------------- | --------- | ----- | ------------- | --- | --- | --- | --- | --- | --- |
| A.2 The       | singular  | value | decomposition |     |     |     |     |     |     |
One of the central tools in tensor-network algorithms is the singular value decomposition
| (SVD). | Any matrix | M of | dimension | D   | A ×D | B can be | written | as  |      |
| ------ | ---------- | ---- | --------- | --- | ---- | -------- | ------- | --- | ---- |
|        |            |      |           |     | M =  | USV†,    |         |     | (17) |
where U has orthonormal columns (U†U = 1), V† has orthonormal rows (V†V = 1),
and S is a diagonal matrix whose entries are non-negative singular values, ordered as
| λ 1 ≥ λ | 2 ≥ ··· ≥ 0. |     |     |     |     |     |     |     |     |
| ------- | ------------ | --- | --- | --- | --- | --- | --- | --- | --- |
In tensor-network language, an SVD is the operation that splits one tensor into two.
The two new tensors are connected by a new internal index, and the dimension of this
index is set by the number of nonzero singular values. In this sense, the SVD not only
factorizesanobject,butalsotellsushowlargethenewbondhastobeinordertorepresent
| the original | tensor | exactly. |     |     |     |     |     |     |     |
| ------------ | ------ | -------- | --- | --- | --- | --- | --- | --- | --- |
The ordering of the singular values is what makes the SVD especially useful. Large
singularvaluescorrespondtothemostimportantcomponentsacrossthechosenbipartition,
while small singular values carry less weight. This gives a controlled way to compress a
tensor: by discarding the smallest singular values, we reduce the dimension of the new
bond while keeping the dominant correlations. This is the basic idea behind truncation.
Wewillusethismechanismrepeatedly. First,itallowsustobuildmatrixproductstates
by splitting a many-body wavefunction site by site. Later, in Section A.6, it provides the
practical rule for reducing bond dimensions after operations that make them grow too
large. The important point is that the SVD gives a systematic way to decide what to keep
| and what   | can be  | safely discarded. |     |           |     |     |     |     |     |
| ---------- | ------- | ----------------- | --- | --------- | --- | --- | --- | --- | --- |
| A.3 Matrix | product | states            | and | operators |     |     |     |     |     |
A matrix product state (MPS) is obtained by applying the SVD idea repeatedly along a
one-dimensional chain. A completely general state on N sites can be written as
X
|     |     |     | |Ψ⟩ | =   | C       | |i  | ,...,i | ⟩,  | (18) |
| --- | --- | --- | --- | --- | ------- | --- | ------ | --- | ---- |
|     |     |     |     |     | i1...iN | 1   | N      |     |      |
i1,...,iN
where the coefficient tensor C contains dN entries. A direct representation is then impos-
sible for large systems. But by sweeping an SVD through the chain, this large tensor can
| be factorized | exactly | into      | a product | of  | smaller | tensors, |      |         |      |
| ------------- | ------- | --------- | --------- | --- | ------- | -------- | ---- | ------- | ---- |
|               |         |           |           | X   | [1]     | [2]      |      | [N]     |      |
|               |         | C i1...iN | =         |     | A       | A        | ···A | ,       | (19) |
|               |         |           |           |     | i1,α1   | i2,α1,α2 |      | iN,αN−1 |      |
α1,...,αN−1
where each tensor carries one physical index i and, except at the boundaries, two virtual
n
| indices | connecting | it to its | neighbours. |     |     |     |     |     |     |
| ------- | ---------- | --------- | ----------- | --- | --- | --- | --- | --- | --- |
The dimensions of the virtual indices are the bond dimensions. They determine how
much correlation the MPS can carry across each cut. If the bond dimension is χ = 1, the
state is simply a product state with no entanglement between the two sides of any cut.
Increasing χ allows the ansatz to represent more entangled states. The factorization itself
is exact if the bond dimensions are allowed to grow as needed; the approximation enters
19

onlywhenweimposeamaximumbonddimensionandtruncatethesmallersingularvalues,
| as discussed | in  | Section | A.6. |     |     |     |     |     |     |
| ------------ | --- | ------- | ---- | --- | --- | --- | --- | --- | --- |
Matrixproductoperators(MPOs)usethesameideabutforoperatorsinsteadofstates.
An MPO is a chain of local tensors with two physical indices, one incoming and one
outgoing, together with virtual indices that connect neighbouring tensors. Applying an
MPO to an MPS generally increases the bond dimension, because the virtual spaces of
the operator and the state combine. This growth is one of the main reasons truncation is
| needed in | practical | tensor-network |     | algorithms.     |     |     |     |     |     |
| --------- | --------- | -------------- | --- | --------------- | --- | --- | --- | --- | --- |
| A.4 Gauge | freedom   |                | and | canonical forms |     |     |     |     |     |
The MPS representation is not unique. On any virtual bond we can insert the identity in
the form 1 = XX−1 and absorb the two factors into the neighbouring tensors,
|     |     |     | A[n] | → A[n]X, | A[n+1] | →   | X−1A[n+1]. |     | (20) |
| --- | --- | --- | ---- | -------- | ------ | --- | ---------- | --- | ---- |
This transformation changes the individual tensors but leaves the physical state |Ψ⟩ un-
changed. This freedom is called gauge freedom. It is extremely useful because it allows us
to choose the tensor representation that is best suited to a given calculation.
The most important choice is the canonical form. By sweeping QR or SVD decomposi-
tions from left to right, tensors can be brought into left-canonical form, meaning that they
satisfy A†A = 1 when contracted over their physical and left virtual indices. Similarly,
a right-to-left sweep gives right-canonical tensors satisfying BB† = 1. If the two sweeps
meet at a chosen site or bond, the MPS is in mixed canonical form: tensors on the left are
left-canonical, tensors on the right are right-canonical, and the center tensor carries the
| nontrivial | weight | of the | state. |     |     |     |     |     |     |
| ---------- | ------ | ------ | ------ | --- | --- | --- | --- | --- | --- |
This form is especially useful for local expectation values. If the orthogonality center is
placedatsiten,thetensorstotheleftandrightcontractwiththeirconjugatestoidentities.
The full network for ⟨Ψ|Oˆ |Ψ⟩ then collapses to a local contraction at the center,
n
|     |     |     |     |       |      | h       | i    |     |      |
| --- | --- | --- | --- | ----- | ---- | ------- | ---- | --- | ---- |
|     |     |     |     | ⟨Ψ|Oˆ |      | Λ[n]†Oˆ | Λ[n] |     |      |
|     |     |     |     | n |Ψ⟩ | = Tr |         | n .  |     | (21) |
Thus, instead of contracting the entire chain, we only need the local tensor carrying the or-
thogonalitycenter. Thisiswhycanonicalformsaresocentralintensor-networkalgorithms:
| they turn   | many | global-looking |     | operations | into    | local   | ones.  |     |     |
| ----------- | ---- | -------------- | --- | ---------- | ------- | ------- | ------ | --- | --- |
| A.5 Schmidt |      | decomposition  |     | and the    | reduced | density | matrix |     |     |
Consider a bipartition cut that divides the chain into two parts: a block A = {1,...,j}
and its complement B = {j +1,...,N}. If we group all physical indices in A into one
composite index and all physical indices in B into another, the wavefunction becomes a
matrix Ψ . Performing an SVD of this matrix gives the Schmidt decomposition,
(A),(B)
χ
|     |     |     | X   |            |     |     | X    |         |      |
| --- | --- | --- | --- | ---------- | --- | --- | ---- | ------- | ---- |
|     |     | |ψ⟩ | =   | s |α⟩ ⊗|α⟩ |     | , s | ≥ 0, | s2 = 1, | (22) |
|     |     |     |     | α A        | B   | α   |      | α       |      |
|     |     |     | α=1 |            |     |     |      | α       |      |
where the states |α⟩ and |α⟩ form orthonormal bases on the two sides of the cut. The
|     |     | A   |     | B   |     |     |     |     |     |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
numbers s are the Schmidt coefficients, and they are exactly the singular values of the
α
| reshaped | wavefunction. |     |     |     |     |     |     |     |     |
| -------- | ------------- | --- | --- | --- | --- | --- | --- | --- | --- |
In an MPS, this structure is already built into the bonds. If the state is brought into
mixed canonical form with the center placed on bond j, the center object is precisely the
20

diagonal matrix of Schmidt coefficients across that cut. The MPS bond dimension χ is
therefore the number of Schmidt values kept on that bond. In this sense, the canonical
| form makes | the entanglement | structure | of every cut | explicit. |     |
| ---------- | ---------------- | --------- | ------------ | --------- | --- |
The reduced density matrix of the block A is diagonal in the Schmidt basis,
X
|     |     | ρ = Tr | |ψ⟩⟨ψ| = s2 | |α⟩ ⟨α| , | (23) |
| --- | --- | ------ | ----------- | --------- | ---- |
|     |     | A      | B           | α A A     |      |
α
so its eigenvalues are simply p = s2. These probabilities determine the entanglement
|     |     | α   | α   |     |     |
| --- | --- | --- | --- | --- | --- |
across the cut, for example through the von Neumann entropy or the Rényi entropies.
| A.6 Optimal | truncation |     |     |     |     |
| ----------- | ---------- | --- | --- | --- | --- |
Truncating a tensor network means reducing some or all of their bond dimensions from χ
|     | χ′  |     |     |     | χ′  |
| --- | --- | --- | --- | --- | --- |
to a smaller value < χ. In the Schmidt basis, the natural way to do this is to keep the
largest Schmidt coefficients and discard the rest. This choice is optimal in a precise sense,
(cid:12) E
given by the Eckart–Young–Mirsky theorem: among all states (cid:12)ψ˜ with Schmidt rank at
(cid:12)
most χ′, the truncated state is the closest one to |ψ⟩ in norm. The corresponding error is
| simply | the discarded | weight, |                                  |      |      |
| ------ | ------------- | ------- | -------------------------------- | ---- | ---- |
|        |               |         | (cid:13) (cid:12) ψ˜ E(cid:13) 2 | X    |      |
|        |               | ε =     | (cid:13)|ψ⟩− (cid:12) (cid:13) = | s2 . | (24) |
|        |               |         | (cid:12) 2                       | α    |      |
α>χ′
This gives a direct and practical criterion for compression: if the discarded Schmidt values
| are small, | the truncation | error is | small. |     |     |
| ---------- | -------------- | -------- | ------ | --- | --- |
For the states of local gapped Hamiltonians, for instance, the Schmidt values typically
decay rapidly, reflecting the limited entanglement expected from the area-law structure
discussedabove. Asaresult,keepingonlythelargestsingularvaluesoftengivesanaccurate
representation with a much smaller bond dimension. In practice, every truncation step
follows this logic. Whenever an operation increases the bond dimension, for example when
applyinganMPOtoanMPS,weperformanSVDandkeeponlytheleadingsingularvalues.
This trades a controlled and usually very small loss of accuracy for a representation that
| remains | computationally    | manageable. |                    |     |     |
| ------- | ------------------ | ----------- | ------------------ | --- | --- |
| A.7 The | reduced transition | matrix      | and its truncation |     |     |
The discussion so far concerned the truncation of a single state |ψ⟩, where the reduced
densitymatrix(RDM)providesthenaturalobjecttodiagonalize. Inthetransversepicture,
the situation is different. The transfer matrix is non-Hermitian, so its dominant left and
rightfixedpointsaredistinct: alefttMPS⟨L|andarighttMPS|R⟩. Theobjectthatmust
be truncated is therefore not an ordinary reduced density matrix, but a reduced transition
| matrix | (RTM) built from | both fixed | points. |     |     |
| ------ | ---------------- | ---------- | ------- | --- | --- |
For a temporal cut that separates the “past” steps ≤ j from the “future” steps > j, the
| RTM is | defined as |     |                  |     |      |
| ------ | ---------- | --- | ---------------- | --- | ---- |
|        |            |     | Tr future |R⟩⟨L| |     |      |
|        |            |     | T =              | .   | (25) |
j
⟨L|R|L|R⟩
Because ⟨L| and |R⟩ are generally not related by Hermitian conjugation, T is itself non-
j
Hermitian. Its eigenvalues can be complex. They still sum to one, but they are not
probabilities, so the usual Schmidt-value interpretation of Section A.5 no longer applies
directly.
This changes the logic of truncation. In the ordinary case, we keep the directions that
bestpreservethenormofasinglestate. Here,thephysicallyrelevantquantitiesareinstead
21

bilinear overlaps of the form ⟨L|O|R|L|O|R⟩. A component of |R⟩ that has little overlap
with ⟨L| does not contribute much to such observables, even if it looks important from
the point of view of |R⟩ alone. Conversely, a component with small weight in |R⟩ can still
| matter | if it has a large | projection | onto ⟨L|. |     |     |     |
| ------ | ----------------- | ---------- | --------- | --- | --- | --- |
Theappropriatetruncationmustthereforeusethetransitionmatrixitself. Themethod
introduced by Carignano, Marimón and Tagliacozzo [CRT24] performs a biorthogonal de-
composition of |R⟩⟨L| and keeps the eigenvalues of largest modulus. These are the di-
rections that dominate the bilinear form and therefore the physical contractions. Since
relevance is measured through the overlap between left and right fixed points, this pro-
cedure can retain fewer states than a standard reduced-density-matrix truncation at the
same accuracy.
The corresponding cost is naturally measured by the generalized Rényi-2 temporal
entropy,
|     |     |     | Sgen(j) | −logTrT2. |     |      |
| --- | --- | --- | ------- | --------- | --- | ---- |
|     |     |     |         | =         |     | (26) |
|     |     |     | 2       |           | j   |      |
This is the entropy we use throughout the thesis to quantify the temporal entanglement
| that controls | the difficulty | of the        | transverse | contraction. |                   |     |
| ------------- | -------------- | ------------- | ---------- | ------------ | ----------------- | --- |
| B Fisher      | Zeros          | and Dynamical |            | Quantum      | Phase Transitions |     |
| review in     | detail this    | section (not  | checked)   |              |                   |     |
This appendix sets out the theoretical background for the dynamical quantum phase
transitions (DQPTs) discussed in the main text. The central idea is an analogy between
non-analyticitiesoftheequilibriumfreeenergyandnon-analyticitiesofadynamicalfreeen-
ergybuiltfromtheLoschmidtamplitude, bothunderstoodthroughthezerosofapartition
| function  | in the complex | plane.             |     |     |     |     |
| --------- | -------------- | ------------------ | --- | --- | --- | --- |
| B.1 Zeros | of the         | partition function |     |     |     |     |
The modern understanding of equilibrium phase transitions as non-analyticities of the free
energy traces back to Lee and Yang [? ? ], who studied the zeros of the grand partition
function as a function of a complexified fugacity. For a finite system the partition function
is a polynomial with strictly complex zeros, so the free energy is analytic on the real axis;
in the thermodynamic limit, however, the zeros condense onto curves, and where such a
curve crosses—or pinches—the real axis, the free energy develops a non-analyticity that
we recognize as a phase transition. Fisher [? ] extended the same picture to zeros in the
complex-temperature plane. The transition is thus encoded not in any single configuration
but in the distribution of zeros and, specifically, in whether they reach the physical axis as
| the system | size diverges. |     |     |     |     |     |
| ---------- | -------------- | --- | --- | --- | --- | --- |
B.2 The boundary partition function and the Loschmidt amplitude
Heyl, Polkovnikov and Kehrein [HPK13] observed that the same machinery applies to
quantum dynamics. Starting from an initial state |ψ 0 ⟩, define the boundary partition
function
|     |     | Z(z) | = ⟨ψ | | e−zH |ψ | ⟩, z ∈ C. | (27) |
| --- | --- | ---- | ------ | ------- | --------- | ---- |
|     |     |      | 0      |         | 0         |      |
For real z = R this is a (thermal-like) expectation value, with |ψ ⟩ acting as a boundary
0
condition on a strip of imaginary-time extent R; for purely imaginary z = it it is the
| Loschmidt | amplitude |      |        |          |            |      |
| --------- | --------- | ---- | ------ | -------- | ---------- | ---- |
|           |           | G(t) | = ⟨ψ | | e−iHt |ψ | ⟩ = Z(it), | (28) |
|           |           |      | 0      |          | 0          |      |
22

the return amplitude to the initial state after a quench. The Loschmidt amplitude is
therefore the analytic continuation of the boundary partition function to the imaginary
axis, and the zeros of Z(z)—the dynamical analogue of Fisher zeros—control its analytic
structure. When, in the thermodynamic limit, these zeros cross the imaginary axis at a
discrete set of points z = it∗, the return amplitude acquires genuine non-analyticities at
n n
the critical times t∗.
n
B.3 Dynamical free energy and dynamical quantum phase transitions
In analogy with the equilibrium free energy density, define the dynamical free energy
1
f(z) = − lim logZ(z), (29)
N→∞ N
and, restricting to real time, the rate function of the return probability
1
(cid:12) (cid:12)2
r(t) = − lim log(cid:12)G(t)(cid:12) = 2Ref(it). (30)
N→∞ N
Adynamicalquantumphasetransitionisapointt∗ atwhichr(t)isnon-analytic—typically
n
a cusp—inherited from a Fisher zero touching the imaginary axis. The name is apt: r(t)
plays the role of a free energy and t the role of a control parameter, with the critical
timesmarkingtheboundariesbetween“dynamicalphases” distinguishedbywhichcomplex
saddle dominates the return amplitude. We emphasize that, like an equilibrium transition,
a DQPT is sharply defined only in the thermodynamic limit; at finite size the zeros sit off
the axis and the cusps are rounded.
B.4 Transfer-matrix realization
The boundary partition function is precisely the object that the transverse contraction
evaluates: written as a contraction in the spatial direction, Z(z) becomes a product of
transfer matrices, and in the thermodynamic limit f(z) = −logµ (z), where µ (z) is the
0 0
leading eigenvalue of the spatial transfer matrix per site. The dynamical free energy is
then analytic wherever µ is the unique dominant eigenvalue and varies smoothly, and
0
it can become non-analytic only where the leading eigenvalue loses that status—that is,
where two eigenvalues cross in modulus and the dominant one is exchanged. This is the
transfer-matrix face of a DQPT: a level crossing in the (complex) spectrum, which in our
modelisorganizedbytheZ symmetryofSection??intoaprotectedcrossingbetweenthe
2
two parity sectors. A comprehensive account of the field is given in the review by Heyl [?
].
23