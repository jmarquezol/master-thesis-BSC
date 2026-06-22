| SciPost | Physics | Codebases      |     |         |     |            |                | Submission |     |
| ------- | ------- | -------------- | --- | ------- | --- | ---------- | -------------- | ---------- | --- |
|         | The     | ITransverse.jl |     | library | for | transverse | tensor network |            |     |
contractions
|     |     |     |     | Stefano | Carignano1⋆, |     |     |     |     |
| --- | --- | --- | --- | ------- | ------------ | --- | --- | --- | --- |
1BarcelonaSupercomputingCenter,08034Barcelona,Spain
⋆stefano.carignano@bsc.es,
Abstract
5202 peS 3  ]hp-tnauq[  1v99630.9052:viXra
Transverse contraction methods are extremely promising tools for the efcient contrac-
tion of tensor networks associated with the time evolution of quantum many-body sys-
tems, allowing in some cases to circumvent the entanglement barrier that would nor-
mallypreventthestudyofquantumdynamicswithclassicalresources. Wepresenthere
the ITransverse.jl package, written in Julia and based on ITensors.jl, containing
severalofthesehigh-levelalgorithms,includingnovelprescriptionsforefcienttrunca-
| tions | of temporal | matrix | product | states. |     |     |     |     |     |
| ----- | ----------- | ------ | ------- | ------- | --- | --- | --- | --- | --- |
Copyrightattributiontoauthors.
ReceivedDate
ThisworkisasubmissiontoSciPostPhysicsCodebases.
AcceptedDate
Licenseinformationtoappearuponpublication.
PublishedDate
Publicationinformationtoappearuponpublication.
Contents
| 1 Introduction |                                          |                                               |            |            |             |         |     |     | 2   |
| -------------- | ---------------------------------------- | --------------------------------------------- | ---------- | ---------- | ----------- | ------- | --- | --- | --- |
| 2 Time         | evolution                                | and                                           | transverse |            | contraction | methods |     |     | 4   |
| 3 Temporal     |                                          | MPOs and                                      | transverse |            | contraction |         |     |     | 7   |
| 3.1            | Modelsandbuildingblocks                  |                                               |            |            |             |         |     |     | 8   |
| 3.2            | tMPOconstructionandrotation              |                                               |            |            |             |         |     |     | 10  |
|                | 3.2.1                                    | TemporalMPOforforwardtimeevolution            |            |            |             |         |     |     | 10  |
|                | 3.2.2                                    | FoldedMPOforexpectationvaluesoflocaloperators |            |            |             |         |     |     | 11  |
| 4 Truncations  |                                          | and                                           | low-level  | algorithms |             |         |     |     | 11  |
| 4.1            | Conventions,overlapsandexpectationvalues |                                               |            |            |             |         |     |     | 12  |
| 4.2            | Canonicalforms                           |                                               |            |            |             |         |     |     | 12  |
| 4.3            | RDMtruncation                            |                                               |            |            |             |         |     |     | 13  |
| 4.4            | RTMtruncation                            |                                               |            |            |             |         |     |     | 14  |
| 5 High-level   |                                          | algorithms                                    | for        | network    | contraction |         |     |     | 16  |
| 5.1            | Powermethod                              |                                               |            |            |             |         |     |     | 16  |
| 5.2            | Lightcone                                |                                               |            |            |             |         |     |     | 19  |
| 5.3            | Caveatsandgeneralconsiderations          |                                               |            |            |             |         |     |     | 20  |
| 5.4            | Somenumericalresults                     |                                               |            |            |             |         |     |     | 21  |
6 Temporal entropies and computational complexity of the transverse contraction 22
1

SciPost Physics Codebases Submission
6.1 Generalizedtemporalentropies 23
6.2 Computationalcomplexity 24
7 Conclusion 25
A Symmetric Singular value and Eigenvalue decompositions 25
A.1 SymmetricSVD 26
A.2 Symmetriceigenvaluedecomposition 26
References 26
1 Introduction
Thestudyofthedynamicsofquantummany-bodysystemsisafascinatingyetextremelychal-
lenging subject in contemporary physics. Several open questions on the topic, such as the
nature of thermalizationor thepossiblebreakdown of ergodicity andthe formation of many-
body localized phases [1,2] that retain memory of their initial conditions, require a detailed
characterizationofthetimeevolutionofquantumsystems,ataskwhich-fromthetheoretical
side- becomes unfeasible using traditional tools due to the exponential cost in representing
wave-functionsformany-particlesystems.
Over the past few years, theoretical condensed matter physics has been revolutionized
by the development of tensor network (TN) algorithms, which have allowed to character-
ize the equilibrium properties of quantum many-body systems with unprecedented accuracy.
Nowadays,methodsbasedontheDensityMatrixRenormalizationGroup(DMRG)[3]formu-
lated in the language of Matrix product states (MPS) are the gold standard for the study of
lower-dimensional quantum many-body systems, and a very active research field has devel-
opedaroundthesetechniques(forreviewsseeeg.[4–6]). Inthiswork,wewillmostlyfocuson
such lower-dimensional systems, particularly 1D spin chains, though many of the techniques
discussedinthefollowingcaninprinciplebeappliedtohigher-dimensionalsystemsaswell.
Traditionally, the expressivity of MPS as a wave-function ansatz for a quantum system is
related to entanglement [7,8], which provides a basis on which an efficient compression of
quantum states can be performed [9]. Since ground states of one-dimensional gapped local
Hamiltonians are known to satisfy an area law for their entanglement, this implies that they
can be efficiently represented using MPS, and optimized by performing local operations with
a computational cost that scales polynomially (instead of exponentially) with the number of
constituentsofthesystem[10]. Since“arealaw”impliesthattheentanglementofabipartition
ofthesystemintwopartsisproportionaltotheareaoftheboundarybetweenthem,thismeans
thatinaone-dimensionalsystem(wheretheareabetweenanytwosegmentsisjusttwopoints)
theentanglementisboundedanddoesnotgrowarbitrarilywiththesizeofthesystem. Avery
different scenario would be given in the case of a “volume law”, for which the entanglement
ofa1Dchainwouldgrowarbitrarilylargeasthenumberofitselementsincreases.
Fromatechnicalpointofview,havingalimitedamountofentanglementresultsinasmall
virtual (often called “bond”) dimension of the tensors forming the MPS, allowing for an effi-
cient compression of the data required to faithfully represent these states. Other well-known
examples of efficient representations using tensor networks are thermal states, which can be
describedbydensityoperators,aswellaslocalHamiltonians(andtheirexponentials),which
allowforacompressedmatrixproductoperator(MPO)form[11–15].
2

SciPost Physics Codebases Submission
Giventhetremendoussuccessofthesemethodsindescribinggroundstateandequilibrium
properties of quantum many-body systems, it is then logical to ask whether they can provide
an efficient representation for studying their time evolution as well. After all, if a system
thermalizes,atypicalscenariosuchasaquenchfromagroundstatetoanotherthermalstate
wouldinterpolatebetweentwolimitscharacterizedbothbylowentanglement.
The answer in this case is less clear: while tremendous progress has been made in devel-
oping highly efficient algorithms based on a Trotterized time evolution such as TEBD [16],
tDMRG[17]andthelike,aswellasvariationalmethodswhichdirectlysolvetheSchroÈdinger
equationonthetangentspaceoftheseMPS[18],ultimatelytheyallrelyonafinalre-compression
of the wave-function which, as before, will depend on the amount of entanglement between
the constituents of the system (see [19] for a comprehensive review). Unfortunately, it is
nowadays well understood that the entanglement entropy of a wave function in a standard
time evolution scenario such as a global quench grows linearly with time [20,21], requiring
in principle an exponential amount of resources to faithfully describe the state in time. This
so-called“entanglementbarrier"constitutesaformidablechallengeforanycurrentTNmethod
basedonconstructingthetime-evolvedwavefunctionforamany-bodysystem.
In recent years, new tensor network algorithms for studying time evolution have never-
theless been proposed. One extremely promising example is given by transverse contraction
methods [22]: the idea is to encode the dynamical evolution of a D-dimensional system in a
D+1-dimensionaltensornetwork,wheretheadditionaldimensionisgivenbytime,andcon-
tract it along the spatial direction using boundary MPS methods (see the following Section).
Improving over the first proposals in this direction, we proposed a novel truncation method
basedonreducedtransitionmatrices,whichallowedtodevelopfurtherintuitiononthecom-
putationalcomplexityofcomputingthedynamicsofquantummany-bodysystemsusingtensor
networks[23]. Thisprescriptionalsoprovidesanaturalconnectionwithconceptsintroduced
in field theory such as generalized temporal entropies [23–25], the physics of open quantum
systems [26] and, importantly, allows to compute expectation values of local operators effi-
ciently[27].
DuringthecourseoftheseinvestigationsheldbyourjointcollaborationatBSC,U.Barcelona
and CSIC1 we began developing a toolkit implementing transverse contraction algorithms,
most of which were lacking a public code implementation. The process culminated with the
creation of the ITransverse.jl package, written in the modern Julia programming lan-
guage [28] and built on top of the ITensors.jl and ITensorMPS.jl libraries [29], which
provide an easy way to write performant code that can be straightforwardly extended to use
GPU acceleration. The idea is to provide a self-contained implementation of most state of
theartalgorithmsrelatedtotransversecontraction,whichallowstheusertoimplementtheir
favorite model and compute in the most efficient way possible time-dependent amplitudes,
expectationvaluesand(generalized)entropies.
InstallationissimplydoneviaJulia’spackagemanager:
julia> using Pkg
julia> Pkg.add(url="https://github.com/starsfordummies/ITransverse.jl.git")
Beforemovingontothedescriptionofthealgorithmsincludedinthelibrary,letusbriefly
reviewtheideasbehindthesetransversecontractionalgorithmsforstudyingthetimeevolution
ofquantummany-bodysystems.
1The people who actively contributed to this line of research are Luca Tagliacozzo, Carlos Ramos MarimoÂn,
AleixBou-Comas,JanT.Schneider,EsperanzaLoÂpezandSergioCerezo.
3

SciPost Physics Codebases Submission
2 Time evolution and transverse contraction methods
The most commonly employed tensor network methods for studying the time evolution of
a quantum many-body system, such as TEBD and TDVP, focus on building the time-evolved
wave-functionofthesysteminacompressedform,correspondingtothefamiliarSchroÈdinger
picture in quantum mechanics. As already anticipated in the previous section, this approach
oftenturnsouttobecomputationallyinefficient,asduringtimeevolutionentanglementtends
togrowwithavolumelaw,renderinganefficientre-compressionusingMPSimpossible.
As is well known, focusing on the time evolution the wave function however is not the
only way to perform time evolution in quantum mechanics. Under many circumstances, we
aremainlyinterestedincomputingexpectationvaluesofa(typicallylocal)operatorO,sayfor
examplethemagnetizationforaspin:
O(t) = ψ(0) U(t) †OU(t)ψ(0) (1)
〈 〉 〈 | | 〉
where ψ(0) denotestheinitialstateand U(t)=exp( iHt)isthetimeevolutionoperator.
Sin|cebot〉htheoperatorandthetime-evolutionoper−atorforalocalHamiltonianH canusu-
allybeefficientlyexpressedinformofmatrixproductoperators(MPOs)[12],anothernatural
optionwouldbetotryandbuildatime-evolvedO(t)=U(t) †OU(t),whichisnothingbutthe
usualHeisenbergpicturefortimeevolutioninquantummechanics. Whilethisapproachturns
out to be computationally advantageous in a few cases, in a generic scenario the “operator
entanglement" which dictates the bond dimension required for an efficient representation of
O(t)alsoturnsouttogrowlinearlywithtime,sothatonceagaintheentanglementbarrierpre-
vents from accessing long-time evolution due to an exponential growth in the computational
resourcesrequired[30,31].
Itmayappearthenthatperforminglong-timeevolutionwithtensornetworksisahopeless
taskduetothegrowthofentanglementbetweenthesitesofthesystem,eitherontheoperator
orthestateside,astimeincreases.
As we will see in the following, tensor networks can however hint towards an additional
way to perform time evolution beyond the traditional SchroÈdinger and Heisenberg pictures.
The key idea for this is to visualize the dynamical evolution of a D-dimensional system as a
D+1-dimensional tensor network, where the additional dimension is given by time. Up to
Trotter errors, this higher-dimensional PEPS of finite bond dimension then encodes the full
dynamicsweareinterestedin;thechallengeisnowtocontractitefficiently.
Wecanstartforexamplebyinspectingthetwo-dimensionalTNassociatedwiththecalcu-
lation of a return amplitude of a product state ψ to itself after some time evolution with a
0
givenHamiltonian,whichisoftenreferredtoas| aL〉oschmidtecho(seeFig.1),
Aψ
0
ψ
0
(t)= ψ(0) U(t)ψ(0) . (2)
〈 | | 〉
The basic ingredients of the network here are the tensors of the initial and final states, rep-
resented as an MPS, and the MPO tensors that make the time evolution operator U(t). For
simplicity,in thefollowing we shallmostlywork withtranslationinvariantsystemsandtime-
independentHamiltonians,sothatalltheseMPOtensorsW inthebulkareidentical.
Now,insteadofevolvingtheinitialproductstatebycontractingrowsoftheU(t)MPOonto
it and building up entanglement among its sites, we can think of contracting this network in
a transverse direction, namely by identifying columns of the TN with states and operators
associatedwithonespatialsiteatdifferenttimesteps[22,32,33]. Thecolumnsattheleftand
right edges of the network will thus now represent a temporal matrix product state (tMPS),
whereas the columns in the middle, which act as spatial transfer matrices (we will usually
refertothemas E,asinFigure1)canbeseenastemporalMPOs(tMPO).Wewillgooverthe
detailsontheconstructionoftheseobjectsinSection3.
4

SciPost Physics Codebases Submission
Figure 1: Two-dimensional (space-time) tensor network representing a time-
dependent return amplitude ψ
0
U(t)ψ
0
: starting from ψ
0
product state, which
we draw at the bottom of the〈ne|twork|(ie〉. time runs upw|ard〉s), we write the time-
evolutionoperatorU(t)astheproductofseveralrowsofMPOs,eachoneassociated
withthetrotterized U(δt)(shadedinblue). Wethencloseagainwiththeconjugate
of the initial state, represented at the top. For a translation-invariant system, the
whole network is given by the repeated product of a single column E, which can be
seenasatransfermatrixthatmovesfromonespatialsitetothenext.
Startingfromtheleftandrightedges,wecanapplythestandardTNmachinerytoperform
thecontractionofthenetworkbyapplyingthetMPOstothetMPSattheedges. Naturally,the
bond dimension of the resulting tMPS can in principle grow exponentially with the number
oftMPO layersappliedtothem,sothat,asusual,someintermediatecompression isrequired
(wewilldiscussthisissueindetailinSection4).
We are then back to the question of whether these temporal states can be efficiently re-
compressed in MPS form at each step. If we turn to the standard DMRG-type algorithms, the
determiningfactorherewillbea“temporal"entanglement,associatedwithagivenelementof
thesystematdifferenttimes[22,32–34]. Whilethephysicalinterpretationofthisquantityis
stillobjectofanactiveinvestigation,fromacomputationalpointofviewthisisawell-defined
question: onejustcomputesthestandardentanglemententropyalongbipartitionsofthetMPS
via its reduced density matrices and checks its scaling with the number of sites, which now
coincides with the number of steps N of time evolution performed: that is, longer chains
t
correspondtolongertimeevolutions.
This prescription was first applied in [22] to compute the expectation value of a local
operator, cf. Equation (1). The tensor network associated with this quantity is shown in
Figure 2 (a): half of the rows are given by the time evolution U(t) MPOs, whereas the other
halfbyitsconjugate U† (t).
Followingintuitionfromsimplertoymodels,inordertoreducetemporalentanglementthe
authorsof[22]suggestedthentofoldthetensornetworkinhalfalongthetemporaldirection,
building composite tensors made of an element of the forward and one of the backwards
evolution MPO (Figure 2 (b)). This construction is reminiscent of the forwards-backwards
Schwinger-Keldysh contour, which is used in real-time field theory calculations. The folding
prescription also allows to build a connection with familiar concepts such as the Feynman-
Vernon influence functional, as pointed out in [34,35]: focusing on one constituent, we can
interpret the left and right vectors built from the contraction to the rest of the network to
its left (or right, respectively) as the description of the “bath", represented by the rest of the
system, at different times. The ability to consider the full spatio-temporal network here can
allowustogoevenfurther,identifyingtemporalMPSasobjectswhichcanbeusedtorepresent
5

SciPost Physics Codebases Submission
Figure2: (a)Two-dimensionaltensornetworkassociatedwiththeexpectationvalue
ofalocaloperator(denotedasaredtensorinthemiddleofthenetwork. (b)Folded
version of the same network: in our convention, the vectorized local operator is at
thetopofthenetwork,whiletheinitialstate(orratherdensitymatrix,inthefolded
picture) is at the bottom. (c) Light cone structure: via the folding operation, all
tensors outside of the causal cone of the local operators reduce to identities, so we
canneglectthemintheactualcomputations.
processtensors,allowingforadirectdescriptionofcorrelationsinspaceandtimeintheunified
languageoftensornetworks[26].
Thefoldedstructureforanexpectationvalueofalocaloperatorallowsalsotoexploitthe
finite speed of propagation of information in the system. By building the TN in a form which
reproducesthelightconeassociatedwiththeoperator,onecaneffectivelyreducethenumber
of tensors required for evaluating its expectation value: all information outside the cone will
not play a role in the dynamics. Specifically, one can see that, upon contraction, all folded
tensorsoutsidetheconecollapsetoidentities(seeFigure2(c)).
HavingdefinedtheboundarytMPSandtMPO,wearenowreadytobuildalgorithmsthat
contract the TN using them. In the ITransverse.jl library we have implemented both a
power method, which can be used to build left and right vectors for systems with an infinite
spatial extension, as well as a light-cone method tailored for the evaluation of expectation
values of local operators, which exploits the causal structure discussed above. We discuss
thesemethodsindetailinSection5.
While all these novel techniques turn out to be advantageous in several situations com-
pared to traditional methods, there are still cases in which temporal entanglement exhibits a
volume law, ie. increases linearly with the time of the evolution considered, effectively pre-
ventingtheiruseforlongtimeevolution.
A closer inspection of the time-evolution TN with a mindset focused on transverse con-
tractionscanneverthelessgiveusevenmoreinformation. Thekeyobservationhereisthat,if
we consider the contraction of the network starting from the sides and construct a tMPS L
associated with the contraction of the left half of the system, as well as a "right" tMPS R 〈for|
the other half, the final result of our calculation will be given by the overlap between|th〉ese
two, L R . The relevant object for our calculation then is not given by the left and right vec-
tors s〈ep|ar〉ately, but rather by their overlap, upon which one can construct a cost function to
optimize when performing a truncated contraction of the network [23]. In turn, this implies
that - rather than thinking about the usual (reduced) density matrices (RDMs) ρ L = L L
andρ R= R R forleftandrightvectorsindividually,oneshouldconsidertransitionma|tr〉ic〈es|
R |L 〉a〈nd| theirreduced(RTMs).
T
As|id〉e〈fr|om providing a possible computational advantage, the method proposed in [23]
allows for a fascinating connection with high energy physics: there the concept of reduced
transition matrices has been recently proposed as starting point to defined generalized en-
6

SciPost Physics Codebases Submission
Figure3: SketchoftherotationimplicitlyperformedinITransverse.jl forindex
labelling: thephysicalindices(i
x
,j x)oftheusualtensorsoftheMPOU(δt)become
virtual indices (α t ,β t) for the temporal MPO, whereas the virtual indices (α x ,β x)
become physical temporal ones (i
t
,j t). The physical Hilbert space of the temporal
degrees of freedom is thus dictated by the virtual dimension of U(δt). Each tMPO
column is made by a stack of rank-4 tensors closed at the top and the bottom (or
rightandleft,afterrotation)bythetrandblstates.
tropies [36–38], which turn out to have a geometrical interpretation in holography [38–42].
Forourcase,therelevantquantityintheoptimizationwouldthenbethegeneralizedtemporal
entropy associated with the left-right contraction of our network. Since the RTMs are not
hermitian by construction, in principle these entropies are complex quantities, so that their
interpretation has to be taken with care. We will discuss the calculation of these entropies,
theirpropertiesandtheirrelationwiththecomputationalcomplexityofcontractingthetensor
networkinSection6.
In the remainder of this work, we will review these aspects associated with transverse
contraction algorithms, exploring along the way the various functionalities implemented in
theITransverse.jl library.
3 Temporal MPOs and transverse contraction
Much like for traditional time-evolution methods such as TEBD, the basic building blocks for
all the algorithms discussed here are the tensors of a matrix product operator, which can be
builtintheusualwaybyfactorizingtheTrotterizedtimeevolutionoperatorU(δt)asaproduct
of local tensors. In order to perform the transverse contraction, we identify columns of these
tensors, defining temporal MPOs (tMPOs) E, as shown in Figure 1. In the same way, at the
boundariesofthenetworkweidentifyleftandrighttemporalMPS,andcontractthenetwork
usingtheseboundarystates. Inthistransversepicture,thelinksofthetimeevolutionMPOthus
becomeassociatedwithaHilbertspace(withphysicaldimensiongivenbythebonddimension
of U(δt)) where temporal degrees of freedom live, while the original physical legs become
linksinthetemporaldirection.
In order to build a connection with the usual notation employed by TN libraries (and to
mostpeople’sintuition),wherewethinkofMPSaschainsextendingfromlefttoright,inthe
following we will often work in a“rotated" picture, namely we considerthe networkas tilted
clockwise by 90 degrees. The top side of the tMPS/tMPO then becomes the right edge of the
rotatedchains. WeillustratethisinFigure3. Ofcourse,thisrotationismerelyavisualization
tool(sothatwecan,eg.,talkaboutleftandrightcanonicalforms),whichhasnoeffectonthe
physicsunderneath.
We start in Section 3.1 describing the basic building blocks of the networks we consider.
7

SciPost Physics Codebases Submission
Figure4: (a): Givenatwo-bodygatesuchasthosethattypicallyappearintheTrotter
decompositionofexponentialsoflocalHamiltonians,wecandecomposethemusing
standardtensorfactorizations. (b): ThetypicalTrotterstructure for U(δt)madeby
the product of even-odd terms can be brought into MPO form after decompositions
and contractions of pairs of the resulting tensors (denoted as shaded regions). As-
suming translation-invariant Hamiltonians, for a open chain of N sites, we end up
with N 2 bulk tensors W sandwiched between a left W and a right W boundary
c l c
tensors.−
Asdescribedintheprevioussection,whendealingwithdynamicsonecantypicallyencounter
two types of scenario: in the first one, we are interested in computing amplitudes of the
from Aφψ(t)= φ U(t)ψ ,sothatwesimplyencodeforwardtimeevolutioninthenetwork,
working with tM〈PS| of N| t 〉 = t/δt sites 2. We discuss the construction of these objects in
Section 3.2.1. In the other case, we are rather interested in computing the expectation value
of a given operator, O(t) ψ = ψ U† (t)OU(t)ψ , and we typically encode both forward
and backwards evolut〈ion in〉 our t〈MP|O representa|tio〉n. The most efficient way of doing this is
usuallytoworkinafoldedrepresentation,thatis,buildingtensorsmadeoftheproductofone
tensorfromtheforwardtimeevolutionoperator,togetherwiththecorrespondingonecoming
fromthebackwardsevolution. WedescribethisinmoredetailinSection3.2.2.
3.1 Models and building blocks
Given a model Hamiltonian H, the building blocks of the corresponding temporal MPOs are
of course the same as the tensors in the MPO representation of the time evolution operator
U(t) = exp( iHt). For a Hamiltonian with local interaction, the idea is typically to perform
a Trotter exp−ansion introducing small time-steps δt, breaking up U(t) = U(δt) and ex-
pressing each U(δt) as the product of local operators such as two-body gates, giving rise to

the familiar brick-wall structure. One can then further decompose these gates and construct
local tensors W, as can be seen in Figure 4, ending up with MPOs for U(δt) built from rows
of local tensors . There are of course several prescriptions for building these tensors, see
eg.[12,14,15],someofwhichareimplementedinITransverse.jl.
ThefactorizationofU(t)providedbytheTrotterexpansionallowstobuildadiscretetem-
poral lattice, where each site is associated with a given instant in time t = n t δt. After the
clockwiserotationdescribedabove,wecanviewthevirtual(bond)linksoftheMPOasphys-
ical temporal sites, with a Hilbert space of a local dimension given by the bond dimension of
U(δt). At the same time, the spatial physical dimension of the initial chain translates into a
virtualdimensionforthetMPOswearebuilding.
The user is free to specify the tMPO tensors for whichever additional model they are in-
terested in3 and use the contraction algorithms described in Section 5. For convenience, in
2Forsimplicity,wecanthinkofcontractingalreadytheinitialandfinalstatesinthefirstandlasttMPOtensors,
whichusuallymakeslifeeasierifwestart/finishwithproductstates.
3Onatechnicalnote,forthepurposeoftransversecontractionandtruncationonreducedtransitionmatrices,
where the hermiticity of the operators involved is not guaranteed, it is often helpful to work at least with sym-
metricobjects,sothatwheneverpossiblewealwaysprefertousesymmetricleft-righttensors(seediscussionsin
Section5.1andSection6.1).
8

| SciPost Physics |     | Codebases |     |     |     |     |     |     |     |     | Submission |     |
| --------------- | --- | --------- | --- | --- | --- | --- | --- | --- | --- | --- | ---------- | --- |
ITransverse.jl we also provide helper functions to initialize the tMPO tensors given the
required parameters for a few commonly used models. Under the hood, they call the appro-
priate build_expH() function for building the tensors of the U(δt) for the model we are
interestedin. Inparticular,weprovideexplicitrepresentationsof U(δt)for
• TheIsingmodel(withtransverseandparallelfields),withHamiltonian
|     |     |     |     |           |     | Jσi | σi +1 | i +hσi |     |     |     |     |
| --- | --- | --- | --- | --------- | --- | --- | ----- | ------ | --- | --- | --- | --- |
|     |     |     | H   | Ising(g)= |     |     | +gσ   |        |     | ,   |     | (3) |
|     |     |     |     |           |     | x   | x     | z      | x   |     |     |     |
−
i
|       |      |                |     |     |   |     |     |     |    |     |     |     |
| ----- | ---- | -------------- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| withσ | andσ | Paulimatrices. |     |     |     |     |     |     |     |     |     |     |
x z
• Thethree-statePottsmodel,definedas
|     |     |     |           |     |     | σ†    | +σ† |      | i+τ† |     |     |     |
| --- | --- | --- | --------- | --- | --- | ----- | --- | ---- | ---- | --- | --- | --- |
|     |     | H   | Potts(g)= |     | J σ |       | σ   | +g(τ |      |     | ,   | (4) |
|     |     |     |           |     |     | i i+1 | i   | i+1  |      | i ) |     |     |
−
i
|          |          |     |     |       |       |     |        |    |     |        |        |       |
| -------- | -------- | --- | --- | ------- | ------ | --- | ------ | --- | --- | ------- | ------ | ----- |
| with the | matrices |     | σ   |         | ωs s s | , ω | ei2π/3 | and | τ   |         | s s+1, | where |
|          |          |     | =   | s=0,1,2 |        | =   |        |     | =   | s=0,1,2 |        |       |
|          |          |     |     |         | | 〉〈   | |   |        |     |     |         | | 〉〈 | |       |
theadditionismodulo3.
|     |     |     |    |     |     |     |     |     |    |     |     |     |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
• TheXXZmodel,whichwecanwriteintermsofthespinoperatorsS (sothatourimple-
i
mentationworksbothforspin1and1/2)
L
|     |     |     |     |     |     | xS x | y   | y    | zS z   |     |     |     |
| --- | --- | --- | --- | --- | --- | ---- | --- | ---- | ------ | --- | --- | --- |
|     |     |     | H   | =   | J S | +S   | S   | +∆S  |        | .   |     | (5) |
|     |     |     | XXZ |     | i   | i +1 | i   | i +1 | i i +1 |     |     |     |
−
i=1
|     |     |     |     |     |   |     |     |     |    |     |     |     |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
AparticularlyusefulformfortheMPOtensorsof U(δt)oftheIsingmodel,whichissym-
metricbothinthephysicalandthevirtuallegs,isgivenbytheprescriptionin[12]andimple-
mentedinthefunction
build_expH_ising_murg(sites::Vector{<:Index}, mp::IsingParams, dt::Number)
wherethestructIsingParamscontainsthecouplingconstantsforthemodel. Itsbulktensors
| associatedwiththetwo-bodytermsfor |     |          |     |          | J =1aregivenby |     |     |           |               |     |     |     |
| --------------------------------- | --- | -------- | --- | -------- | -------------- | --- | --- | --------- | ------------- | --- | --- | --- |
|                                   |     |          |     | c o s(   | δ t )         |     |     | isin( δ t | ) c o s (δt)σ |     |     |     |
|                                   | W   | b u lk   |     |          |                |     |     |           |               | x   |     | (6) |
|                                   |     | =        |     |          |                |     | −   |           |               |     |     |     |
|                                   |     | I s in g |     | isin ( δ | t )c o s(δt)σ  |     |     | s in      | ( δ t )       |     |     |     |
|                                   |     |          |    |          |                | x   |    |           |              |    |     |     |
|                                   |     |          | −   |          |                |     |     | −         |               |     |     |     |
could
While in principle we use the same prescription for the Potts model, the resulting
tensors do not have left-right symmetry. To get it, we can instead use a symmetric SVD de-
composition[43,44] of thetwo-bodygates(thesesymmetricdecompositionsare particularly
useful,sowerecallthemexplicitlyinAppendixA.1),ensuringthatweendupwithsymmetric
tensors. ThisimplementationofthePottsmodelisimplementedinthefunction
build_expH_potts_symm_svd(sites, mp::PottsParams, dt::Number)
In addition to these models, for more generic Hamiltonians, we also provide a generic
MPOformforU(t)uptosecond-orderTrotter,followingtherecipein[15],whichcanbebuilt
startingfromtheformofthelocalMPOtensorsofagivenHamiltonian4. Ingeneralthisform
isnotsymmetric,butcanbestillusedinourcode.
4ThankstoJanT.Schneiderforsharingtheimplementation.
9

SciPost Physics Codebases Submission
3.2 tMPO construction and rotation
InITransverse.jl wedirectlyprovidehelperfunctionstoinitializethetMPOtensorsgiven
the required parameters, for the models discussed above. Since we usually rely on Trotter-
ized versions of exponentials of Hamiltonians, we require as input a time-step dt, a function
expH_func()tobuildtheindividualtensors(seeprevioussection)andastructcontainingthe
model parameters for the Hamiltonian. When providing these quantities to our constructors,
we build the local tensors W for both the real time evolution operator U(δt) as well as its
imaginary time counterpart, U( iδt). We can specify if we want to incorporate a few steps
of imaginary time in our evolut−ion via the nbeta parameter. Finally, since often in our cal-
culations we consider translation-invariant initial states that do not change throughout our
simulations,weprovide the possibilityto specify oneof these repeating tensors fortheinitial
state,whichinourconventionsendsupatthebottomendofthetMPO(or,afterthe90degrees
rotation described above, at the left edge - hence the “bottom-left / bl " name). All together,
theparametersforthetMPOconstructionarecontainedinthe
struct tMPOParams
dt::Number
expH_func::Function
mp::ModelParams
nbeta::Int
bl::ITensor # bottom -> left(rotated)
3.2.1 Temporal MPO for forward time evolution
Having specified the input parameters, we can use the function fw_tMPO to build the tMPO
associated with the "forward” time evolution U(t), with the boundaries contracted with the
initialandfinal(spatial)states. Forasystemwhichistranslationallyinvariantinspace,thisis
theonlyingredientnecessarytobuildthecontractionofthefullnetwork.
Since we may want to deal with finite systems (or at least use the boundaries of a finite
system as initial guess for the power method described later), under the hood we start by
constructingtheMPOforU(δt)forafinitesystem,andstoringthethreetensorswhichwillbe
used as building blocks for all our network, namely the left and right edges of U(δt) and the
bulktensor,weshallrefertothemasW,W andW respectively,seeFigure4. Forconvenience,
l r c
we also compute the tensors of U( iδt) associated with imaginary time evolution: these
can be used to regularize by cooling−down initial states before performing the evolution by
replacing the first few W of the tMPO with the corresponding W . These tensors are stored
im
in
struct FwtMPOBlocks
Wl::ITensor
Wc::ITensor
Wr::ITensor
Wl_im::ITensor
Wc_im::ITensor
Wr_im::ITensor
tp::tMPOParams
rot_inds::Dict
andtherot_indsdictionaryprovidesthemappingbetweenthelegsoftheunrotatedandthe
rotated W tensors, associating the virtual links of the spatial MPO to the physical sites of the
rotatedtMPOandviceversa.
10

SciPost Physics Codebases Submission
Figure 5: Construction of the folded tMPO tensors, made of pairs W and W†. To
build the WW we take their outer product and merge the legs in pairs, resulting
in a single tensor which has squared dimensions compared to the original one. For
convenience,wedrawheresomearrowsinordertoshowtheorderingofthelegsin
theproduct.
3.2.2 Folded MPO for expectation values of local operators
As we already mentioned, when dealing with expectation values of local operators it is often
convenienttoworkwithfoldedtMPOs,whosebuildingblocksWW W W†aresimplygiven
bytheouterproductofaW tensorfromtheforwardtimeevolution≡cont⊗ourtogetherwiththe
correspondingonefromthebackwardsevolutionpath(seeFigure5).
Thesetensorsareconstructedandstoredinthe
struct FoldtMPOBlocks
WWl::ITensor
WWc::ITensor
WWr::ITensor
WWl_im::ITensor
WWc_im::ITensor
WWr_im::ITensor
rho0::ITensor
tp::tMPOParams
rot_inds::Dict
whichhasthesamestructureasFwtMPOBlocks,thoughherewealsoallowtheusertostore
the(folded)tensorrho0associatedwiththeinitialstate,whichcanbeseenasthelocaltensor
partofadensitymatrix.
TheseblocksarebuiltandputtogetherinafoldedtMPObythefunction
folded_tMPO(b::FoldtMPOBlocks, ts::Vector{<:Index}; fold_op, outputlevel)
where fold_op, which defaults to an identity, specifies the operator we want to put at the
foldingpointofthetMPObetweenforwardandbackwardstimeevolution,whileoutputlevel
determineswhetheradditionaldebugginginformationsshouldbeprinted.
4 Truncations and low-level algorithms
Before moving to the high-level algorithms performing transverse contraction, we start by
introducingthetruncationschemeswhichareemployedunderneath.
As customary, most of the transverse contraction algorithms rely on a basic operation,
namely the application of tMPOs to left and right tMPS. In turn, a truncation of the bond
dimension of the resulting tMPS is required in order to keep the computational cost under
control,andweimplementinourlibraryseveralprescriptionstailoredforthiskindofproblem.
Afterbrieflyreviewing ourconventions in Section 4.1andcanonicalformsinSection 4.2,we
11

SciPost Physics Codebases Submission
Figure6: (a): MPStensorsA incanonicalformreducetoidentitieswhencontracted
i
with their conjugate (in the appropriate direction, depending on the form of the
isometry, so we can talk about left and right canonical tensors). (b): We can think
of a “generalized” canonical form involving two distinct MPS by imposing that the
contraction of their tensors A ,B associated with the same physical site reduces to
i i
an identity in a similar way. Here and in the following we use the convention that
shapesofthesameformandinlight/darkershadesofthesamecolordenoteatensor
anditscomplexconjugate,respectively.
discussthetruncationbasedonreduceddensitymatricesinSection4.3andtransitionmatrices
inSection4.4.
4.1 Conventions, overlaps and expectation values
Thefullcontractionofthe2Dtensornetworkassociatedwiththetimeevolutionprovidesthe
dynamicalquantitiesweareinterestedin. Ifwethinkofcontractingthelefthalfofthesystem
into a single “left tMPS”, and the right half into a “right tMPS”, the full contraction is simply
givenbytheirproduct.
When thinking about states, we are used to computing overlaps of two “kets” represent-
ing vectors ψ , φ via the scalar product φ ψ , which involves the complex conjugation
(denoted he|re〉wi|th〉a ) going from the ket t〈o t|he〉bra. For example, if we take a qubit state
∗
φ = (a,b), we will have φ = (a ,b ) T. In our transverse contraction setup, the quantity
∗ ∗
|we〉are interested instead is〈si|mply given by the product of the two halves of the system, no
conjugation needs to be made. We take this into account by considering the left half of the
systemdirectlyasa“bra” L ,whereastherighthalfisreferredtoas R 5. Theiroverlap L R
doesnotimplythusanyco〈n|jugation. Forthis,weprovideinITrans|ve〉rse.jl thefun〈cti|on〉
overlap_noconj(L::MPS,R::MPS)
which,given L and R ,computes L R withoutanyconjugation6.
Another s〈cen|ario|w〉hich can be〈of|in〉terest is the contraction of the left and right halves
of the system into tMPS with the exception of a central tMPO column E . As an example, it
O
could be that we were able to compute L and R and wanted to use them to compute the
expectationvalueofdifferentlocalopera〈to|rsont|he〉centralsiteofthesystem,whichcouldbe
included in different tMPOs (hence the label “O", as in operator). The full contraction of the
networkcanbewrittenthenas L E R ,whereagainnoconjugationismade. Weprovidethe
O
helperfunctionexpval_LR(L,〈Eo|,R|)〉tocomputethisoverlapinanefficientway.
4.2 Canonical forms
Canonicalformsarewellknowntomosttensornetworkpractitionersasoneofthemostpow-
erfulinstrumentsallowingtoperformlocaltruncationswhicharegloballyoptimal,aswellas
to greatlysimplify calculations andimprove thenumerical stability of algorithms (seeeg.[4]
5Ifthenetworkwassymmetricleft-right,wewouldthensaythat L = R ∗ .
6This is of course equivalent to inner(dag(L),R) using ITen〈so|rs’〈in|ner() function, but we do it in a
moreefficientwaywhichnaturallyavoidsthedoubleconjugation.
12

SciPost Physics Codebases Submission
Figure7: Thecommonlyemployedtruncationbasedonreduceddensitymatrices: if
we are in the appropriate mixed canonical form, we can access the spectrum of the
RDMatagivencutbysimplyperformingasingularvaluedecompositionofthetensor
in the orthogonality center, which is equivalent to the eigenvalue decomposition of
a similar matrix to the RDM. This in turn allows to perform the optimal truncation
(in norm) of the full state by performing a series of local operations, since all the
relevant information is contained in the tensor in the orthogonality center. In the
case shown here, the orthogonality center is in the central (square) tensor, and the
wholeRDMspectrumcanbeobtainedbyitsSVD.
for a review). Given a MPS, we can bring it to a mixed canonical form on a given site by
ensuringthatallitstensorstotheleftofachosenorthogonalitycentersiteareleftisometries,
whilethosetotherightarerightisometries. ThiscanbeseenastheconditionthateachMPS
tensor contracted with its conjugate in the appropriate direction reduces to a simple identity
(seeFigure6). Inthisform,wecanthinkoftherelevantinformationonthefullstateasbeing
contained in the orthogonality center tensor, while all the others provide a unitary change of
basis from the physical Hilbert space of the full chain into the virtual legs of the MPS. If the
orthogonalitycentercorrespondstothefirstsite,theMPSissaidtobeinrightcanonicalform,
whereasitisleft-canonicalifthecenterisonthelastsite.
For later convenience, we can define here a less conventional condition on our tensors,
this time involving two MPS sharing the same physical sites, which in our case will end up
being the tMPS L and R . We impose namely that the product of a tensor of the first MPS
with the corresp〈on|ding o|n〉e of the second, taken again in the appropriate left-right direction,
reduces to an identity [23]. In fact, we can even consider the symmetric case R = L , in
∗
whichcasethistranslatesintoaconditiononasingleMPS.Withanabuseofna|m〉ing,〈we|call
tMPSsatisfyingthisconditionasbeinginageneralizedcanonicalform7
4.3 RDM truncation
The simplest algorithm that can be employed for compressing L and R is familiar to each
TN practitioner, and involves truncating over the largest eigenv〈ec|tors o|f t〉he reduced density
matricesbuilt fromeachtMPS,takenseparately. Afterputting the L and R into the appro-
priate mixed canonical forms, this is equivalent to performing sing〈ula|r valu|e〉decompositions
(SVDs) on the tMPS tensors, see Figure 7. This truncation, as we all know, provides the best
approximationin2-normforthe L and R vectorsindividually. Ofcourse,thisisastandard
truncation used in most MPS libr〈ari|es, an|d〉we can simply call ITensor’s apply() function to
performit.
Onemayarguehoweverthatthistypeoftruncationisnotthebestfortheproblemweare
considering: for non-hermitian objects such as the spatial transfer matrices considered here,
theresultoftheoptimizationturnsouttobegauge-dependent. Toseeit,wecanconsiderthe
2DTNmadebyaseriesoftMPOcolumnsEsandwichedbetweentheleftandrightvectors,and
applyatransformation(whichcanbechosenforsimplicitytobelocal)viathetensorsX,which
leavesthenetworkinvariant,redefining E aswellas L and R ,seeFigure8. While E issent
into a similar matrix, if we try to build its left and rig〈ht|dom|in〉ant vectors using their RDM ρ
7Notethatthisformwillbeill-definedifthetwovectorsareorthogonal.
13

SciPost Physics Codebases Submission
Figure 8: Effects of gauge transformations: we depict in (a) in a schematic way
the 2D network we want to contract, represented as the contraction of a left vector
with a right one (optionally one can include transfer matrices E in the middle and
consider that the network contraction will be given by the overlap of the dominant
eigenvectorsofE). Forsimplicityandwithoutlossofgeneralityweconsiderhereonly
twotime-steps. Wecanalwaysperformhereagaugetransformationbyinsertingan
identityintheform XX 1,with X anyinvertiblematrix,redefiningtheleftandright
−
temporal MPS. (b) While in principle this has no effect on the network contraction,
thespectraoftheRDMfor L and R willdependonthisgaugechoice. (c)Wecan
seeinsteadthatthespectru〈m|ofthe|re〉ducedtransitionmatrixbuiltfrom L and R
isgauge-independent. 〈 | | 〉
Figure 9: We can compute efficiently the singular values of the RTM if we are in
the appropriate canonical form (we refer to Figure 6 for a recap of our graphical
notations). As an example, here we start in left canonical form, and show explicitly
how-thankstothepropertiesofthecanonicaltensors-wecancomputetheSVDby
simplycontractingthelasttensorsoftheleftandrighttMPS.
toperformourtruncation,weseethattheapplicationofX hasastrongeffectontheirspectra.
Crucially,X willingeneralnotbeunitary,andinfactcanevenbeseenassomeimaginarytime
evolution which can reduce the entanglement of the boundary states, possibly even turning
them into product states. This reduction of rank of the RDM does not translate however into
amoreefficientrepresentationoftheTN:thetwostatesmaybecomeincreasinglyorthogonal
astherelevantcontributionstotheiroverlapgetshiftedthetailsofsingularvalueswhichget
truncatedintherequiredcompressionprocedure,leadingtoaninevitable(anduncontrolled)
lossofprecisioninthefinalresult(see[45]foramoredetaileddiscussionofthisissue).
4.4 RTM truncation
InordertocircumventthelimitationsoftheRDMtruncationdescribedabove,itcanbeuseful
to take a step back and recall that ultimately we are interested in the full contraction of the
2D TN associated with the dynamical quantity we are looking at. This change of perspective
suggests that another prescription for truncating the left and right tMPS is possible: one can
namely try to find the best approximation for the overlap L R , rather than focusing on the
twovectorsindividually. Afirstproposalinthisdirectionwa〈sg|iv〉enby[33],wheretheauthors
proposedanalgorithm8foratranslation-invariantsymmetricsystembasedontruncatingover
thelargestsingularvaluesofthematrix R R (orequivalently L L ).
∗ ∗
| 〉〈 | | 〉〈 |
8Wesometimesrefertoitasthe“Hastings”truncation,afterthefirstauthors’name.
14

SciPost Physics Codebases Submission
Figure 10: Example of truncation algorithm based on reduced transition matrices,
startingfromthelastsitesofthetwoinputtMPS(werefertoFigure6forarecapof
our graphical notations). We first bring the two tMPS into left canonical form, then
perform a right sweep, building environments which at each cut will have the same
singular values as the corresponding RTM. At each step of the sweep, we truncate
by inserting the appropriate isometries and redefining the tensors of the two tMPS,
keepingonlythelargestsingularvalues. Aswemovealongthetemporalchains,we
diagonalizethepreviousenvironments(possiblyevenreducingthemtoidentitiesby
appropriatelymultiplyinganddividingbythesquarerootsofthesingularvalues,so
effectively bringing the two tMPS in “generalized” canonical form. Here we show
explicitly the first step of this truncation with the tensor operations involved, the
procedureisthenrepeatedforallsites.
Inspired by that recipe and realizing its connection with the concept of transition matri-
ces, in [23] we developed a truncation algorithm based on reduced transition matrices built
from the left and right vectors encoding the full network (including the operators, if we are
interestedincomputingexpectationvalues):
R L
t
=Tr
t| 〉〈 |
, (7)
T L R
〈 | 〉
whicharenon-hermitianmatriceswithunittrace9. ThefinalcontractionofthefullTNwillbe
given by the overlap of the L and R obtained by contracting over the left and right halves
ofthenetwork. 〈 | | 〉
WhileideallywewouldliketotruncateovertheeigenvaluesoftheseRTMs,inanalogyto
what one does with reduced density matrices, non-hermitian matrices will not have a series
of positive definite (or even real) eigenvalues, so that a truncation prescription here is not
obvious. Instead, in order to provide the best low-rank approximation for our objects we can
choose to to truncate on the largest singular values of these RTM, a prescription which is
alwaysrobustevenfornon-hermitianobjects.
The algorithm, which is sketched in Figure 10, goes as follows [23]: we start by bringing
both L and R individuallyintocanonicalformintheusualway(weprovidetheoptiontodo
itfrom〈 |leftto|ri〉ghtandviceversa,whichcanleadtoinequivalentresultsifthenetworkisnot
top-downsymmetric),letusassumeinthefollowingthatwebringitintoleftcanonicalfrom.
Wethenbeginourtruncation: westartconstructingtherightenvironment bycontracting
N
E
together the last tensors of L and R , call them L[N] and R[N], where N is the number of
temporal sites of the two tM〈PS|. We|th〉en proceed to compute the SVD of this environment,
which–thankstothecanonicalform–correspondstotheSVDofthereducedtransitionmatrix,
as canbe seen in Figure 9 (where we employed the usualtrick that similar matrices have the
same eigenvalues, eig(UAU 1 )=eig(A)). We now use the resulting isometries to diagonalize
−
the environment, transforming the tensors L[N] and R[N]. In fact, we can even transform
9Inthisframework,thealgorithmproposedin[33]canbeseenasatruncationbasedonRTMinthesymmetric
case L = R ,withouttheinsertionofanoperator.
〈 | | 〉
15

SciPost Physics Codebases Submission
into an identity by multiplying L[N] and R[N] each by the square root of the inverse of
N
E
the singular values: in doing so, we bring the left and right tMPS in a form such that the
contractionsofthe L[i]andR[i]totherightofourworkingsiteresultinidentities,seeFigure6
(b). WethenproceedalongthetwotMPS,buildingenvironmentsanddiagonalizingthemfor
eachtemporalsiteandtruncatingovertheirlargestsingularvalues.
We can ask ourselves at this point why this particular prescription arises only in the case
of the transverse contraction, or, equivalently, to what it would correspond in the case of the
usualtimeevolutionaÁlaSchroÈdingerinthetemporaldirection. Uponinspectingthenetwork
structure,wecanseethatinthatcasethe"bottom"and"top"vectorswhichwouldreplace R
and L intheconstructionofthetransitionmatricesarepreciselywhatwillendupbecomi|ng〉
ψ(t〈 ) |and its conjugate ψ(t). So in that case the transition matrix would have exactly the
|same〉shape as the densit〈y mat|rix for ψ(t) , so that the two prescriptions (RTM and RDM)
coincide. | 〉
Bothtypesoftruncationaretypicallycontrolledbythetypicalparametersoneencounters
in a TN calculation: the cutoff below which we discard singular values, and the maximum
bonddimensionwecanafford. Forconvenience,wepassthemaroundusingthe
struct TruncParams
cutoff::Float64
maxbondim::Int64
direction::String # "left" or "right"
wherefortheRTMmethodwealsoallowtheusertospecifythedirectionofthetruncation
sweep,whichcanleadtoinequivalentresultforasymmetricnetworkssuchasthoseconsidered
here(seethediscussioninSection5.3).
5 High-level algorithms for network contraction
InITransverse.jl weprovidetwohigh-levelalgorithmsforperformingthetransversecon-
traction of the 2D network associated with the dynamics of our system. The first is a generic
powermethodwithdeterminestheleftandrightdominanteigenvectorsofthetransfermatrix
by repeatedly applying it to starting tMPS. As such, this method works best for a homoge-
neous system, and it can be applied both for the construction of time-dependent amplitudes
as well as the calculation of expectation values. We describe it in Section 5.1. The second
methodisthelightconealgorithmproposedin[46],whichisbestsuitedforthecalculationof
expectation values of local operators in the thermodynamic limit, as we show in Section 5.2.
AfterdiscussingsomecaveatsandgoodpracticesinSection5.3,weshowcasesomenumerical
resultsobtainedusingthesemethodsinSection5.4.
5.1 Power method
Thepowermethodisastandardalgorithminlinearalgebrawhichallowstocomputeinmost
cases the dominant eigenvectors of a given matrix, namely those with the largest absolute
value (assuming that a gap is present, otherwise the method will be bound to the subspace
of degenerate dominant eigenvectors). The idea is simply to start from a trial vector and
repeatedlyapplythematrixweareinterestedindiagonalizingtoit,theresultwillincreasingly
converge to the dominant eigenvector we are interested in. For non-normal matrices, the
procedurewillneedtobedoneseparatelyfortheleftandrightvectors.
Considering now MPS and MPO as vectors and matrices which we can efficiently multi-
ply, this method provides a way to effectively determine left and right eigenvectors if we are
interestedindeterminingthedynamicsofatranslation-invariantsystem. ThisboundaryMPS
16

SciPost Physics Codebases Submission
Figure11: Powermethodalgorithm: foratranslation-invariantinfinitenetwork,the
full network contraction will be given by the overlap of the left and right dominant
vectors of one column of the network, which acts as a transfer matrix E. (a): By
repeatedlyapplyingthistMPOtoaninitialtMPS φ ,weendupprojectingoverits
0
dominant vectors, with a convergence rate dictat|ed b〉y its gap. (b): At each step of
thepowermethodweapplythetMPO E totheleftandrighttMPS,andre-compress
them,repeatinguntilconvergenceisreached(seemaintextfordetails).
methodiscommonlyemployedinthecontractionof2Dnetworks,andwepresenthereourim-
plementation tailored for the contraction associated with the dynamics. In practice, for most
cases,aslongasweevolveforafinitenumberoftime-stepsthetransfermatricesweconsider
haveagap[24],guaranteeingarelativelyquickconvergenceofthemethod.
LetusnowbrieflydescribetheactualimplementationsprovidedinITransverse.jl. The
powermethodparametersarecontrolledviathe
struct PMParams
truncp::TruncParams
itermax::Int64
eps_converged::Float64
increase_chi::Bool
opt_method::String
The power method runs for a maximum of itermax iterations, or stops earlier if conver-
gence has been reached. The latter is estimated as the variation of the entropy ∆S, which is
computed via RDM or RTM, depending on the optimization method used (as determined by
opt_method, see below and Section 4) between the tMPS computed in subsequent iterations
of the power method: if ∆S < eps_converged, we say that the power method has con-
verged10. In some cases, the bond dimension reached in the intermediate steps of the power
methodissignificantlylargerthantheonerequiredforthefinalconvergedtMPS.Infact,ifwe
haveagoodguessonthebonddimensionrequiredforthefinalresult,itissometimespossible
to specify a maxbondim close to that for the algorithm, which in this case will typically need
moreiterationstoconverge,buteachoneatalowercomputationalcost. Inanycase,wecan
try to slow down the growth of this bond dimension in the hope that convergence is reached
beforelargerdimensionsarehit: thisisgovernedbytheincrease_chiparameter,which-if
settotrue-limitsthemaximumbonddimensionintheinitialiterationsofthepowermethod,
whilelettingitgraduallyincreasetowardstheend.
We provide different functions to perform the power method, depending on the problem
athand. Thefirstoptionisprovidedbythefunction
powermethod_both(in_mps::MPS, in_mpo_L::MPO, in_mpo_R::MPO,
pm_params::PMParams; flip_R::Bool=false)
10Otheroptions,suchasconsideringthefidelitybetweenthetMPSinastepandthatofthefollowingoneare
possible,ofcourse.
17

| SciPost Physics | Codebases |     |     |     |     | Submission |     |
| --------------- | --------- | --- | --- | --- | --- | ---------- | --- |
which, starting from an initial guess init_mps on both sides, at each step applies a column
in_mpo_L to the left and a in_mpo_R to the right, then performs the truncation either indi-
viduallyonthenew L and R usingRDM,oroptimizingtheRTM R L ,asdescribedabove.
Sobothleftandrigh〈tv|ector|sa〉reoptimized,hencethename. | 〉〈 |
Whenever possible, working with systems that have a left-right symmetry with respect to
the center allows to greatly simplify the contraction procedure, since we only need to update
oneofthetwotMPSandobtaintheotherbysimpletransposition. Weimplementthesymmet-
ricpowermethodforthiscaseinthefunction
powermethod_sym(in_mps::MPS, in_mpo::MPO, pm_params::PMParams)
whichtakesasinputasingleMPSin_mps,appliesthein_mpotoitandtruncatestheresulting
MPS, repeating until convergence is reached, as before. We can specify in the
pm_params
opt_method = RDM if we want to use the usual truncation based on the reduced density
matrix, or opt_method = RTM, in which case, since the problem is symmetric, we compute
the SVD of the symmetric environments in a Autonne-Takagi form (see Appendix A.1). In
fact, in this symmetric case one can also compute efficiently the (complex!) eigenvalues of
the RTM (we discuss this in more detail in Section 6.1), so we also provide an experimental
truncation which discards the eigenvalues with the smallest modulus, which can be selected
| viaopt_method | = RTM_EIG. |     |     |     |     |     |     |
| ------------- | ---------- | --- | --- | --- | --- | --- | --- |
Suppose now that we are interested in the contraction of a network associated with the
expectationvalueofalocaloperatorinaninfinitesystem. Inthetransversepicture,wecansee
the final result as the contraction L ENE EN R , for N so that the initial guesses
|     |     |     | N O | N   |     |     |     |
| --- | --- | --- | --- | --- | --- | --- | --- |
1 1
fortheboundaries L and R 〈ar−en|otimport|ant.〉Here E→isatMPOcolumnwithoutany
|     | N   | N   |     | 1   |     |     |     |
| --- | --- | --- | --- | --- | --- | --- | --- |
d〈en−ot|esthe|col〉umn(s)containingadditionsoflocaloperatorsbetweenthe
| operator,while | E   |     |     |     |     |     |     |
| -------------- | --- | --- | --- | --- | --- | --- | --- |
O
forwardsandbackwardscontour. Forthiskindofsetup,wealwaysworkinthefoldedpicture,
whichprovidesthemostefficientrepresentationofourproblem.
Here one might be tempted to optimize left and right vectors without the central column
E , and include it only at the end to evaluate the expectation value. This however will not
O
work: withoutanoperator,thecontractionofthefoldedtensorsreducestoaseriesidentities,
infactthewholenetworkwouldjustreducetothetrivialoverlap ψ U† (t)U(t)ψ ψ ψ .
|               |            |            |                |             | 0   | 0 =                   | 0 0 |
| ------------- | ---------- | ---------- | -------------- | ----------- | --- | --------------------- | --- |
|               |            |            |                | sta〈tes|for |     | and| R〉, lo〈sing|all〉 |     |
| The RTM-based | truncation | would then | return trivial | product     |     | L                     |     |
|               |            |            |                |             |     | 〈 | | 〉               |     |
theinformationonthedynamics.
Instead,onecanusetheprescriptionimplementedinthefunction
function powermethod_op(in_mps::MPS, in_mpo_1::MPO, in_mpo_O::MPO,
pm_params::PMParams)
which goes as follows: for each step, starting from the current L and R we build an up-
|               | byapplyingin_mpo_1toit,andarightone〈R| |            |        |     | by|ap〉plyingin_mpo_O |             |     |
| ------------- | -------------------------------------- | ---------- | ------ | --- | -------------------- | ----------- | --- |
| datedlefttMPS | L                                      |            |        |     |                      |             |     |
|               | 1                                      |            |        |     | O                    |             |     |
|               | Th〈e|resulting                         |            |        |     | t|hen〉used           |             |     |
| to the right. |                                        | transition | matrix | R L | is                   | to optimize | L , |
|               |                                        |            | T      | O 1 |                      |             | 1   |
which is taken as the new L , while R is discar∼de|d. T〉h〈e s|ame procedure is then repe〈ated|
O
the other way around: we〈bu|ild L | an〉d R , optimize their RTM and obtain the updated
|     |     |     | O 1            |     |     |     |     |
| --- | --- | --- | -------------- | --- | --- | --- | --- |
|     |     | W〈  | t|hen re|pe〉at |     |     |     |     |
R , which will be the new R . e until convergence is reached. This method,
1
w| hi〉ch works for a generic ca|se〉 without any left-right symmetry, is employed if we pass the
pm_params.opt_method=RTM_LRparameter.
IfthetMPOsaresymmetricleft-right,weshouldbeabletoperformtheupdateonlyonone
bytransposition11.
ofthetwovectors,say L ,andobtainimmediately R Thisprescriptionis
| usedifwepassthepm_〈pa|rams.opt_method=RTM|_R〉 |     |     |     | parameter. |     |     |     |
| --------------------------------------------- | --- | --- | --- | ---------- | --- | --- | --- |
Finally, for convenience we also allow here to perform the power method by using the
common truncation based on the RDM of R at each iteration. In this case, which is used
| 〉
11Notethatthisis*not*thesameasperformingasymmetricupdateoftheform RR ,herewestillupdatethe
| overlap LO1R |     |     |     |     | 〈 | | 〉   |     |
| ------------ | --- | --- | --- | --- | --- | --- | --- |
| 〈 |          | 〉   |     |     |     |     |     |     |
18

SciPost Physics Codebases Submission
Figure 12: After transverse contraction, the left and right halves of the light cone
TN can be represented as tMPS, which can be built iteratively for each timestep, as
discussedinthemaintext.
Figure13: UpdatingtherighttMPSusingthelightconealgorithm. Ateachstep,we
apply to the current tMPS a tMPO which has one more temporal site on top (built
fromthetensorW ),andobtainanewtMPSwithoneadditionalsiteandlargerbond
r
dimension. WethenproceedtotruncatethetMPS.Inthisway,thetMPSforN time-
t
stepscanbeobtainedbytheonefor N 1stepswithjusttheapplicationofasingle
t
tMPO. −
whenpm_params.opt_method=RDM,thein_mpo_Oinputisunused,andthepowermethod
isperformedusingin_mpo_1only.
Since the tMPOs are usually not unitary operators, at each step of the power method we
normalize back the tMPS, in order to avoid losing precision over many iterations. The most
consistent way to do it is to enforce that the overlap L R = 1, but in practice normalizing
individually L L = RR =1 works as well. Another〈p|os〉sibility is to normalize the overlap
LO 1R befo〈re|tr〉unc〈ati|ng〉 for Rnew .
〈 | 〉 | 〉
5.2 Light cone
Initially proposed in [46] (see also [47]), the light cone method provides an efficient way of
buildingtheleftandrighttMPSassociatedwiththeexpectationvalueofalocaloperator. The
ideabehindthisiterativemethodispreciselytoexploitthecausal(or"light")coneassociated
withsuchoperatortoworkwithasignificantlyreducednumberoftensor[48],seeFigure12.
OnecanstartwithaninitialsetofleftandrighttMPS L and R fortime t ,whichcanbe
1 1 1
given by a single tensor if we start from t 1 = δt, and〈fro|m the|re〉build the tMPS associated
withthetimeevolutionforallsubsequenttimes. Inpractice,inordertobuildthetMPS L
i+1
and R ,oneappliestothematMPOwithanadditionaltimesite,whicheffectivelyex〈tends|
i+1
them| (Fig〉ure13):
R i R i+1 = E i+1 R i , L i L i+1 = L i E i+1 , (8)
| 〉 → | 〉 | 〉 〈 | → 〈 | 〈 |
then a truncation is performed, in a similar way to the power method. The advantage of the
light cone is that, given the tMPS describing the folded left and right vectors at a time t ,
1
we are able to build from them the tMPS at a later time t > t constructively from them by
2 1
applyinglayersoftMPO,withoutneedingtore-contractthefullnetworkuntilconvergence.
19

SciPost Physics Codebases Submission
Given the causal structure dictated by the local operator, the expectation value obtained
will match the value in the thermodynamic limit. In its simplest version, for each additional
timestep the algorithm extends the system by one spatial site, allowing for a description of
physical velocities up to v max = 1/δt. In fact, for most cases this prescription is excessive,
as physical Lieb-Robinson velocities usually fall below this value. One can then work with
narrowerlightcones,extendingtheconeinspaceonlyeveryfewtimesteps[46].
InITransverse.jlwehaveimplementedthefunction
init_cone(tp::tMPOParams, n::Int)
whichinitializesthetMPSforntimestepsandagivensetoftMPOParams,withoutperforming
anytruncation. Thefunction
run_cone(psi::MPS, b::FoldtMPOBlocks, cp::ConeParams, nT_final::Int)
thenevolvestheinputtMPSpsiuntilnT_finaltimesteps,withthetimeevolutionMPObuilt
fromthefoldedblockscontainedinthestructb.
Theparametersofthealgorithmsarecontainedinthe
struct ConeParams
truncp::TruncParams
opt_method::String
optimize_op::Vector
which_evs::Vector{String}
which_ents::Vector{String}
checkpoint::Int
vwidth::Int
where optimize_op allows to define the operator we want to optimize for, using the same
strategyasinthefunctionpowermethod_op()describedabove,ie. truncatingusingRTMscon-
tainingtheoperator(infact,theparameteropt_methodisthesameasforthepowermethod).
Sincethealgorithmevolvesintimestepbystep,wecancomputetheexpectationvalueoflocal
operators and (generalized) temporal entropies at each iteration to follow their time evolu-
tion. Weallowtheusertospecifywhichquantitiestocomputeviatheparameterswhich_evs
andwhich_ents. Finally,weallowtoaccountforphysicalvelocitiessmallerthan1/δt viathe
vwidth parameter (defaulting to 1), which provides the number of timesteps after which the
coneshouldbeextendedbyonespatialsite.
5.3 Caveats and general considerations
We can see the algorithms presented here as falling into two groups: the first one is tailored
forsystemswithasimplestructure,liketheLoschmidtechoforatranslationinvariantsystem,
where the whole network is basically described by a single column which is repeated indefi-
nitely. In that case, the strategy for performing updates is straightforward: we simply apply
oneofthesecolumnstotheleftandtheright,andthensubsequentlytruncate. Asweargued
intheabovesections,thebesttruncationhereshouldbebasedontransitionmatrices.
The second group is tailored to problems which, in a way or another, break translation
invariance,suchastheexpectationvaluesoflocaloperators. Takeasanexamplethelightcone
for an operator , which is built step by step by inserting longer tMPO columns in the center,
“pushing" the previous tMPS outwards by one site, that is, we are growing the system from
the center outwards. A central observation made in [23] is that, if we begin our truncation
sweepfromthesideoftheoperatoratagivenstep,wecandirectlyconnectwiththeoperator
space entanglement of the time-evolved O. It turns out however that this truncation is not a
good basis for the light cone algorithm described above: the issue is that for each time-step,
we are optimizing the tMPS associated with the expectation values of the local operator at
20

SciPost Physics Codebases Submission
Figure14: Comparisonbetweentraditional(TEBD)andtransversecontractionmeth-
ods (both RDM and RTM-based, see main text) for an infinite Ising chain. (a): Ex-
pectation value of σ at mid-chain for an integrable Ising model, with a maximum
z
bonddimensionof64: weseethatTEBDstartsdeviatingfromtheexactresultalready
aroundt 8J,whiletransversecontractionmethodsfaithfullyreproducethecorrect
value. (b)≈: BonddimensionrequiredbythealgorithmsforanintegrableIsingmodel
with a fixed cutoff 10 8. All transverse contraction algorithms exhibit a logarithmic
−
scaling ofthebond dimension,hinting atefficientsimulability ofthedynamics with
tMPS, compared to the exponential scaling of TEBD. (c): Expectation value of σ
z
for non-integrable Ising with transverse and parallel fields, enforcing a maximum
bonddimensionχ =128. Againweseethattransversecontractionmethods,partic-
ularly those based on RTM truncation, are able to capture the correct results with a
limitedχ,whereasTEBDfailsalreadyaround t =5J. (d): Bonddimensionfornon-
integrableIsing. HereRDM-basedtransversetruncationalsoexhibitsanexponential
growth,whiletheRTMtruncationseemstoshowamorefavorablescaling.
the given time we are considering. There is however no guarantee that the optimized tMPS
whichwillgivethecorrectexpectationvalueatatime t canbeusedasgoodstartingpointfor
building the tMPS for a system at t > t. In other words, we might be truncating over a non-
′
optimal quantity, discarding information which will be required at later steps. A workaround
to this is to truncate starting from the side of the initial state by setting direction="left"
in TruncParams. In the absence of an operator and for a symmetric case, this boils down to
thetruncationproposedin[33]andemployedalreadyin[46].
5.4 Some numerical results
Whilewedeferamorecompletediscussiononthecomputationalcomplexityofthetransverse
contractiontothenextsection,asanappetizerletusshowcaseheresomenumericalresultsof
thealgorithmsfor the simple caseof theexpectation valueof alocaloperator,built using the
lightcone. WewillcomparebothRDMandRTMtruncationswiththeresultsobtainedwitha
standardTEBDalgorithm,takenasprototypeofthetraditionaltimeevolutionmethodswhich
evolvethewave-functionintime.
21

SciPost Physics Codebases Submission
In Figure 14 we show results for both the integrable transverse field Ising model, as well
as the non-integrable version with parallel field (Equation (3)). For the latter, we choose a
parametersetwhichshoulddescribeachaoticdynamics,presentingaconsiderablechallenge
forthestudyoftime-evolution.
We start by showing the expectation value of the σ operator at mid-chain, compared to
z
the exact value12 in the various algorithms for a fixed maximum bond dimension, showing
thattransversemethodsprovideamuchsuperiorperformance.
We also compare the bond dimensions required by the algorithms if we impose a cutoff
of 10 8 on the truncated singular values. The difference is striking in the integrable case
−
(b), where we see that the bond dimension χ required by the light cone algorithm with RTM
truncationexhibitsalogarithmicgrowth,remainingbelowχ =100evenbeyondt =10J. The
RDM truncation also seems to have a logarithmic behavior, albeit with a larger prefactor, in
contrastwiththeexponentialgrowthofregularTEBD(notethelogscale). Thenon-integrable
case (d) looks qualitatively different: the bond dimension for transverse methods based on
RDMtruncationscalesexponentially,whereastheRTM-basedmethodalsogrowssignificantly
butseemssub-exponential.
Many other numerical results obtained with ITransverse.jl can be found in the liter-
ature [23–27]. We encourage interested users to check the examples/ folder in the package
repository,whereweprovideseveralscriptstoplayaroundwithdifferentmodelsandtrunca-
tions.
6 Temporal entropies and computational complexity of the trans-
verse contraction
After presenting all the tools for transverse contraction, it is now time to discuss when these
algorithms can be more useful compared to the tried and true SchroÈdinger and Heisenberg
evolutions. Whileingeneralthequestionisstillopen,therearealreadyseveralresultswhich
providesomeinsightinthisdirection,andwebrieflyreviewtheminthissection.
We need not stress here the relevance of the concept of entanglement entropy and all its
implications. Fromthepointofviewofatensornetworkpractitioner,oneofthemostimpor-
tant relationships is the one between the entropy of a state and bond dimension required for
its faithful description: the higher the entanglement, the larger the bond dimension. In turn,
this means that we can only represent efficiently states with small entanglement: for one-
dimensional chains we usually think about states fulfilling an area law (where the entangle-
ment,andthusthebonddimension,remainbounded),orwithatmostlogarithmicviolations
toit,correspondingtoapolynomialresourcecostintermsofmemoryforthetensors.
Inthissense,theentanglemententropyofastate,builtfromitsreduceddensitymatrices,
isaclearindicatorofthecomputationalcomplexityofrepresentingfaithfullythatstateusing
MPS. A natural construction for our case is then to compute the entropy of left and right
temporal states, and expect that such a temporal entropy will provide an estimate of the cost
of computing the dynamical properties we are interested in - after all, the final result will be
givenbytheoverlap L R .
We have seen alre〈ad|y〉however in Section 4.3 that for our case one has to be careful, and
even if the temporal states taken individually can be deceivingly simple, their overlap may
be not accurately reproduced. We can argue instead that the proper object reflecting the
complexity of representing the dynamics is encoded in the transition matrices defined above,
whicharebuiltfromthetwovectorstogether[23].
12Forthenon-integrablecasewetreatas"exact"thevalueobtainedwithTEBDusingthelargestbonddimension
wewereabletouse,χ=4096inourcase.
22

SciPost Physics Codebases Submission
Figure 15: We can compute the generalized Renyi-2 entropy efficiently by simply
constructingtwocopiesoftheRTM(L1,R1andL2,R2,respectively)andcontracting
them appropriately. Here we draw schematically the TN contraction for Tr 2 for
Tt
acutinthemiddleofthetemporalchains.
 
Let us then introduce more formally the concept of generalized temporal entropies, and
showandhowcantheybecalculatedusingITransverse.jl.
6.1 Generalized temporal entropies
The idea of generalized entropies has been proposed in recent high-energy physics literature
[36–38,49] discussing transition matrices between different states. In our case, these will be
theleftandrighttemporalstates,sothatwecantalkaboutgeneralizedtemporalentropies[38,
39,42]. Ifwenowconsider,forexample,thegeneralizedVonNeumannentropy,
gen
S
1
(t)= Tr
T t
log
T t
,
−
we usually think about this expression in terms of the eigenvalues of the RTM . Assum-
t
T
ing that this non-hermitian matrix canbe diagonalized, these eigenvalues will in principle be
complex,leadinginturntocomplex-valuesentropies.
Fromacomputationalpointofview,thereisanadditionalcomplicationhere: theRTMare
exponentially large matrices in the number of time-steps of the time evolution, so that exact
diagonalizationisobviouslyunfeasibleforthembeyondveryshorttimes. Forthecalculationof
reduceddensitymatrices,theMPSmachineryallowsforanefficientdiagonalizationthanksto
the canonical form, which greatly simplifies the problem: the singular values extracted from
an MPS tensor in the orthogonality center are precisely the square roots of the eigenvalues
of the RDM at that cut. Unfortunately, we cannot use straightforwardly the same trick for a
RTM, built from two different vectors: even after a gauge transformation in a “generalized"
canonical form, the tensors on the two sides will be different, so that we cannot rely on a
similaritytransformationtosimplifyourproblemhere.
Oneexceptiontothisisgivenonceagainbythesymmetriccase,inwhich L = R . There,
wecandiagonalizetheRTMwith(complex)orthogonalmatrices[50](seeA〈pp|end|ix〉A.2),so
that we are able to build an environment which is similar (ie. has the same eigenvalues) to
theRTMwewanttodiagonalize.
Thissymmetricdiagonalizationisimplementedinthefunction
diagonalize_rtm_symmetric(psiL::MPS;
bring_left_gen::Bool=true, normalize_eigs::Bool=true,
sort_by_largest::Bool=true, cutoff::Float64=1e-12)
Even if we cannot diagonalize RTMs directly in the non-symmetric case, we can still de-
fine higher-order entropies, which can be computed efficiently using MPS. The first obvious
examplesarethe(generalized)purityandtheα=2Renyientropy,Figure15
S 2= logTr
Tt
2 . (9)
−
 
23

SciPost Physics Codebases Submission
Theevaluationofthesegeneralizedentropiesisprovidedbythefunctions
gen_tsallis2(psi::MPS, phi::MPS)
gen_renyi2(psi::MPS, phi::MPS)
which, given two MPS psi and phi, build the appropriate contraction at each temporal cut
andcomputeefficientlytherequiredtraces.
6.2 Computational complexity
Now that we introduced the definition of generalized temporal entropies, we can finally ask
thequestiononwhatisthecomputationalcomplexityofthetransversecontractionalgorithms
describedhere.
There are several works in the literature studying the scaling of the “standard" temporal
entanglement which comes into play in the RDM truncation, particularly in the context of
discrete(Floquet)dynamics,seeeg.[22,46,51–55].
Let us discuss here instead the scaling of generalized temporal entropies, which - as we
argued-shouldbemorecloselyrelatedtothe dynamicalproblemsweareinterested inhere,
suchastheevaluationoftime-dependentexpectationvaluesandamplitudes.
The RTM truncation described in Section 4.4 allows for a first connection to determine
the computational cost of computing expectation values of operators. By performing this
truncation procedure, one can see that the rank of the RTM involved in the contraction is
upperbounded[23]bytherankofthematrixencodingtheoperatorspaceentanglementen-
tropy[56]. Intheworstcasescenario,thisisstillaquantitythatcouldexhibitalineargrowth
withtime[57,58],ie. avolumelawcorrespondingtoanexponentialentanglementbarrier.
Inordertogetabetterintuitiononthecomputationalcomplexityoftheproblemandthe
role of generalized temporal entropies, it is useful to take a step back and focus on a simpler
TN, namely the one associated with a time-dependent amplitude such as a Loschmidt echo
(Figure1). In[24],thiskindofsetupforaquenchofaninfinitesystemevolvedwithacritical
Hamiltonian was mapped to a path integral on a strip, allowing to exploit the machinery of
conformal field theory (CFT) to predict the behavior of generalized temporal entropies. As
longasananalyticalcontinuationfromeuclideantorealtimeevolutioncanbeperformed,the
growth of these entropies in this context dictates the compressibility of the relevant left and
right temporal MPS, much in the same way as what happens for the spatial entanglement of
critical ground states. Furthermore, CFT predicts a logarithmic growth of generalized entan-
glementatcriticality,hintingatapolynomialcostandthusefficientsimulabilityofamplitudes
usingtMPS.
Recently,anumericalexplorationbeyondthisanalyticalregimehasshownthatthegener-
alizedtemporalentropiesstillexhibitalogarithmicgrowthaftercrossingdynamicalquantum
phase transitions - that is, times at which physical quantities exhibit discontinuities [59,60].
Hinting at some universal properties for large times, this suggests that there is no change in
complexity related to these transitions [27], and that transverse contraction methods might
provideaccesstothelong-timedynamicsofquantumsystemsinanefficientway.
Finally,werecallthat,strictlyspeaking,forageneralnon-symmetriccasetheRTMtrunca-
tionisbasedonasingularvaluedecomposition. Therankoftheobjectsinvolvedinthatcase
couldthenberelatedtotheconceptofSVDentropies,whichhavealsobeenrecentlydiscussed
in the literature [61,62]. A detailedstudy of thesequantities would alsobe of extremeinter-
est,andcouldpossiblyprovidemorerigorousboundsforthecomputationalcomplexityofour
methods.
24

SciPost Physics Codebases Submission
7 Conclusion
Transverse contraction methods provide an efficient way to study the dynamics of quantum
many-bodysystemswithtensornetworks,allowing(atleastinsomecases)tocircumventthe
exponential entanglement barrier associated with the conventional SchroÈdinger and Heisen-
berg time evolutions. By considering the D+1-dimensional tensor network associated with
thetemporalevolutionofaD-dimensionalsystemandintroducingboundary“temporal"states,
thenetworkcontractioncanbeperformedinthespatialdirection. Thecomplexityofthecon-
traction for these boundary MPS methods is related to the (generalized) temporal entangle-
ment, built from reduced transition matrices. For one-dimensional quantum systems, there
isevidencethatmanytemporalstatesassociatedwithcomputingtime-dependentamplitudes
satisfyanarealawforthisgeneralizedentanglement,oratmostalogarithmicgrowth,hinting
atthepossibilityofanefficientrepresentationusingtemporalmatrixproductstates. Thiseffi-
cient representation however does not come completely for free: the typical objects entering
into transverse contraction algorithms for time evolution are complex and non-hermitian, so
thatadditionalcaremustbetakenforthenumericalstabilityandreliabilityofthealgorithms.
We presented here the ITransverse.jl julia package, which provides state of the art
algorithmsbaseduponthetransversecontraction forthedynamicsofone-dimensionalquan-
tum many-body systems. The package provides both high-level routines to perform the full
dynamicalcalculationofexpectationvaluesandamplitudes,includingapowermethodanda
light cone algorithm, as well as several functions to efficiently perform truncations based on
reducedtransitionmatricesandcomputegeneralizedentanglements. Allthesefunctionalities
makeitapowerfultoolbox,whichallowstostudydifferentnumericalproblemsinanefficient
wayusingthesenovelalgorithms.
Whilethealgorithmspresentedherearemeantforthetimeevolutionofone-dimensional
systems, many of the ideas reviewed here can be extended to higher-dimensional problems -
infact,someinitialworksinthisdirectionhavealreadybeenpresentedinthisdirection[63].
WelookforwardtoimplementingthesemethodsinITransverse.jl inthenearfuture.
Acknowledgements
ThealgorithmsincludedinITransverse.jl havebeendevelopedtoinvestigateallsortsof
physical questions related with the dynamics of quantum many-body systems. Most of this
work has been performed in close collaboration with Luca Tagliacozzo, together with Aleix
Bou-Comas, Jan T. Schneider, Sergio Cerezo RoquebruÂn, Esperanza LoÂpez, Jacopo de Nardis,
Guglielmo Lami and Carlos Ramos MarimoÂn. Special thanks to Luca Tagliacozzo, Aleix Bou-
Comas,JanT.SchneiderandMariCarmenBanÄulsforhelpfuldiscussionsandfeedbackonthe
manuscript.
Funding information SC acknowledges his AI4S fellowship within the “GeneracioÂn D” ini-
tiativebyRed.es,MinisterioparalaTransformacioÂnDigitalydelaFuncioÂnPuÂblica,fortalent
attraction(C005/24-EDCV1),fundedbyNextGenerationEUthroughPRTR.
A Symmetric Singular value and Eigenvalue decompositions
The transition matrices which appear in transverse contractions are typically complex and
non-hermitian operators, requiring extra care in their treatment. Luckily, in some cases we
can consider models with a left-right symmetry in the corresponding MPO tensors - in this
25

SciPost Physics Codebases Submission
case, these transfer matrices are symmetric, ie. E = ET (where ET denotes only transposi-
tion, without complex conjugation), a property which allows us to perform a few very useful
decompositions,asdescribedbelow.
A.1 Symmetric SVD
A complex symmetric matrix allows for a special singular value decomposition reflecting this
symmetry,oftenreferredtoasanAutonne-Takagidecomposition[43,44]. Givenasymmetric
matrix M, we can start by performing the standard SVD M = USV†, truncating as usual over
thesmallestsingularvaluesifnecessary. Fromthis,weconstructthe(block-diagonal,ordiag-
onal if there are no degenerate singular values) matrices Z = U†V = V†U = ZT, and from
∗ ∗
theseand U
Z
=UZ. Thisallowsustoreachthefollowingdecompositionusinganisometry
anditstranspose
M =U
Z
SU
Z
T, (A.1)
A.2 Symmetric eigenvalue decomposition
Symmetric matrices also allow for an eigenvalue decomposition using (complex) orthogonal
matrices[50],namely
M OΛOT, OT =O 1. (A.2)
−
≈
Togetthisform,wefirstobtaintherighteigenvectorsV
R
ofM,MV R=V
R
Λ,fromwhichwecan
buildthetransformationmatrixO byfirstconstructingtheblock-diagonalmatrix G=(V
R
TV R)
andthencomputingO=V
R
G
−
1/2.
This decomposition is particularly useful as it defines a similarity transformation which
allows us to diagonalize the reduced transition matrices in a symmetric case using low-rank
representations,asdiscussedinSection6.1.
References
[1] P.Sierant,M.Lewenstein,A.Scardicchio,L.VidmarandJ.Zakrzewski, Many-bodylocal-
ization in the age of classical computing*, Reports on Progress in Physics 88(2), 026502
(2025), doi:10.1088/1361-6633/ad9756.
[2] R. Nandkishore and D. A. Huse, Many-body localization and thermalization in quan-
tumstatisticalmechanics, AnnualReviewofCondensedMatterPhysics6(1),15(2015),
doi:10.1146/annurev-conmatphys-031214-014726.
[3] S. R. White, Density matrix formulation for quantum renormalization groups, Phys. Rev.
Lett.69,2863(1992), doi:10.1103/PhysRevLett.69.2863.
[4] U. SchollwoÈck, The density-matrix renormalization group in the age of matrix product
states, Ann.Phys.(N.Y).326(1),96(2011), doi:10.1016/j.aop.2010.09.012.
[5] S.-J. Ran, E. Tirrito, C. Peng, X. Chen, L. Tagliacozzo, G. Su and M. Lewenstein, Ten-
sor Network Contractions: Methods and Applications to Quantum Many-Body Systems,
Springer International Publishing, ISBN 9783030344894, doi:10.1007/978-3-030-
34489-4(2020).
[6] M. C. BanÄuls, Tensor Network Algorithms: A Route Map, Annual Review of Condensed
Matter Physics 14(1), 173 (2023), doi:10.1146/annurev-conmatphys-040721-022705,
2205.10345.
26

SciPost Physics Codebases Submission
[7] L. Amico, R. Fazio, A. Osterloh and V. Vedral, Entanglement in many-body systems, Re-
viewsofModernPhysics80(2),517(2008), doi:10.1103/RevModPhys.80.517.
[8] N.Laflorencie, Quantumentanglementincondensedmattersystems, PhysicsReports646,
1(2016), doi:10.1016/j.physrep.2016.06.008.
[9] J.I.Cirac,D.Perez-Garcia,N.SchuchandF.Verstraete,Matrixproductstatesandprojected
entangled pair states: Concepts, symmetries, theorems, Rev. Mod. Phys. 93(4), 045003
(2021), doi:10.1103/RevModPhys.93.045003, 2011.12127.
[10] J.Eisert,M.CramerandM.B.Plenio,Colloquium: Arealawsfortheentanglemententropy,
ReviewsofModernPhysics82(1),277–306(2010), doi:10.1103/revmodphys.82.277.
[11] F. Verstraete, J. J. GarcÂõa-Ripoll and J. I. Cirac, Matrix Product Density Operators: Sim-
ulation of Finite-Temperature and Dissipative Systems, Phys. Rev. Lett. 93(20), 207204
(2004), doi:10.1103/PhysRevLett.93.207204.
[12] B. Pirvu, V. Murg, J. I. Cirac and F. Verstraete, Matrix product operator representations,
NewJournalofPhysics12(2),025012(2010), doi:10.1088/1367-2630/12/2/025012.
[13] C. Hubig, I. P. McCulloch and U. SchollwoÈck, Generic construction of efficient matrix
productoperators, Phys.Rev.B95,035129(2017), doi:10.1103/PhysRevB.95.035129.
[14] M. P. Zaletel, R. S. K. Mong, C. Karrasch, J. E. Moore and F. Pollmann, Time-evolving
a matrix product state with long-ranged interactions, Phys. Rev. B 91, 165112 (2015),
doi:10.1103/PhysRevB.91.165112.
[15] M. V. Damme, J. Haegeman, I. McCulloch and L. Vanderstraeten, Efficient higher-
order matrix product operators for time evolution, SciPost Phys. 17, 135 (2024),
doi:10.21468/SciPostPhys.17.5.135.
[16] G.Vidal, EfficientSimulationofOne-DimensionalQuantumMany-BodySystems, Physical
ReviewLetters93(4),040502(2004), doi:10.1103/PhysRevLett.93.040502.
[17] S. R. White and A. E. Feiguin, Real-Time Evolution Using the Density Ma-
trix Renormalization Group, Physical Review Letters 93(7), 076401 (2004),
doi:10.1103/PhysRevLett.93.076401.
[18] J. Haegeman, J. I. Cirac, T. J. Osborne, I. Piˇzorn, H. Verschelde and F. Verstraete, Time-
dependentvariationalprincipleforquantumlattices, Phys.Rev.Lett.107,070601(2011),
doi:10.1103/PhysRevLett.107.070601.
[19] S. Paeckel, T. KoÈhler, A. Swoboda, S. R. Manmana, U. SchollwoÈck and C. Hubig, Time-
evolutionmethodsformatrix-productstates, arXiv:1901.05824[cond-mat,physics:quant-
ph](2019), 1901.05824.
[20] P. Calabrese and J. Cardy, Entanglement entropy and quantum field theory, Journal of
StatisticalMechanics: TheoryandExperiment06,002(2004), doi:DOI:10.1088/1742-
5468/2004/06/P06002;eprintid: arXiv:hep-th/0405152.
[21] M. Fagotti and P. Calabrese, Evolution of entanglement entropy following a quantum
quench: Analytic results for the xy chain in a transverse magnetic field, Phys. Rev. A
78,010306(2008), doi:10.1103/PhysRevA.78.010306.
[22] M. C. BanÄuls, M. B. Hastings, F. Verstraete and J. I. Cirac, Matrix Product States for Dy-
namicalSimulationofInfiniteChains, PhysicalReviewLetters102(24),240603(2009),
doi:10.1103/PhysRevLett.102.240603.
27

SciPost Physics Codebases Submission
[23] S.Carignano,C.R.MarimoÂnandL.Tagliacozzo, Ontemporalentropyandthecomplexity
of computing the expectation value of local operators after a quench, Physical Review
Research6(3),033021(2024), doi:10.1103/PhysRevResearch.6.033021, 2307.11649.
[24] S. Carignano and L. Tagliacozzo, Loschmidt echo, emerging dual unitarity
and scaling of generalized temporal entropies after quenches to the critical point,
doi:10.48550/arXiv.2405.14706(2024),2405.14706.
[25] A. Bou-Comas, C. R. MarimoÂn, J. T. Schneider, S. Carignano and L. Tagliacozzo, Mea-
suring temporal entanglement in experiments as a hallmark for integrability (2024),
2409.05517.
[26] S. Cerezo-RoquebruÂn, A. Bou, J. Schneider, E. LoÂpez, L. Tagliacozzo and S. Carig-
nano, Spatio-temporal tensor-network approaches to out-of-equilibrium dynamics bridg-
ing open and closed systems, Frontiers in Quantum Science and Technology 4 (2025),
doi:10.3389/frqst.2025.1568471.
[27] S. Carignano, G. Lami, J. D. Nardis and L. Tagliacozzo, Overcoming the entanglement
barrierwithsampledtensornetworks(2025),2505.09714.
[28] J.Bezanson,A.Edelman,S.KarpinskiandV.B.Shah,Julia: Afreshapproachtonumerical
computing, SIAMReview59(1),65(2017), doi:10.1137/141000671.
[29] M. Fishman, S. R. White and E. M. Stoudenmire, The ITensor Software Li-
brary for Tensor Network Calculations, SciPost Phys. Codebases p. 4 (2022),
doi:10.21468/SciPostPhysCodeb.4.
[30] T. Prosen and I. Piˇzorn, Operator space entanglement entropy in a transverse Ising chain,
PhysicalReviewA76(3),032316(2007), doi:10.1103/PhysRevA.76.032316.
[31] J. Dubail, Entanglement scaling of operators: A conformal field theory approach, with a
glimpseofsimulabilityoflong-timedynamicsin1+1d, JournalofPhysicsA:Mathematical
andTheoretical50(23),234001(2017), doi:10.1088/1751-8121/aa6f38, 1612.08630.
[32] A. MuÈller-Hermes, J. I. Cirac and M. C. BanÄuls, Tensor network techniques for the com-
putation of dynamical observables in 1D quantum spin systems, New Journal of Physics
14(7),075003(2012), doi:10.1088/1367-2630/14/7/075003, 1204.5080.
[33] M. B. Hastings and R. Mahajan, Connecting Entanglement in Time and Space:
Improving the Folding Algorithm, Physical Review A 91(3), 032306 (2015),
doi:10.1103/PhysRevA.91.032306, 1411.7950.
[34] A.Lerose,M.SonnerandD.A.Abanin, InfluenceMatrixApproachtoMany-BodyFloquet
Dynamics, PhysicalReviewX11(2),021040(2021), doi:10.1103/PhysRevX.11.021040.
[35] M.Sonner,A.LeroseandD.A.Abanin, Influencefunctionalofmany-bodysystems: Tempo-
ralentanglementandmatrix-productstaterepresentation, AnnalsofPhysics435,168677
(2021), doi:10.1016/j.aop.2021.168677.
[36] Y. Nakata, T. Takayanagi, Y. Taki, K. Tamaoka and Z. Wei, Holographic Pseudo Entropy,
PhysicalReview D103(2), 026005 (2021), doi:10.1103/PhysRevD.103.026005, 2005.
13801.
[37] A. Mollabashi, N. Shiba, T. Takayanagi, K. Tamaoka and Z. Wei, Pseudo-Entropy
in Free Quantum Field Theories, Physical Review Letters 126(8), 081601 (2021),
doi:10.1103/PhysRevLett.126.081601.
28

SciPost Physics Codebases Submission
[38] K. Doi, J. Harper, A. Mollabashi, T. Takayanagi and Y. Taki, Pseudo Entropy in dS/CFT
and Time-like Entanglement Entropy, Physical Review Letters 130(3), 031601 (2023),
doi:10.1103/PhysRevLett.130.031601, 2210.09457.
[39] K.Doi,J.Harper,A.Mollabashi,T.TakayanagiandY.Taki,Timelikeentanglemententropy,
doi:10.48550/arXiv.2302.11695(2023),2302.11695.
[40] M. P. Heller, F. Ori and A. Serantes, Geometric interpretation of timelike entanglement
entropy, Phys.Rev.Lett.134,131601(2025), doi:10.1103/PhysRevLett.134.131601.
[41] M.P.Heller,F.OriandA.Serantes, Temporalentanglementfromholographicentanglement
entropy (2025),2507.17847.
[42] T. Takayanagi, Essay: Emergent holographic spacetime from quantum information, Phys.
Rev.Lett.134,240001(2025), doi:10.1103/pg4r-fy8n.
[43] L. Autonne, Sur les matrices hypohermitiennes et sur les matrices unitaire, Annales de
l’UniversiteÂdeLyon,NouvelleSeÂrieI(Fasc.38),1(1915).
[44] T. Takagi, On an algebraic problem related to an analytic theorem of carathÂeodory and
fejÂer and on an allied problem in linear algebra, Japanese Journal of Mathematics 1, 83
(1925), OriginalworkonAutonne-Takagifactorization.
[45] W.Tang,F.VerstraeteandJ.Haegeman, Matrixproductstatefixedpointsofnon-hermitian
transfermatrices,Phys.Rev.B111,035107(2025),doi:10.1103/PhysRevB.111.035107.
[46] M. FrÂõas-PeÂrez and M. C. BanÄuls, Light cone tensor network and time evolution, Physical
ReviewB106(11),115117(2022), doi:10.1103/PhysRevB.106.115117, 2201.08402.
[47] T.EnssandJ.Sirker, Lightconerenormalizationandquantumquenchesinone-dimensional
Hubbard models, New Journal of Physics 14(2), 023008 (2012), doi:10.1088/1367-
2630/14/2/023008, 1104.1643.
[48] M.B.Hastings,Light-conematrixproduct,JournalofMathematicalPhysics50(9)(2009),
doi:10.1063/1.3149556.
[49] S. Murciano, P. Calabrese and R. M. Konik, Generalized entanglement entropies in two-
dimensionalconformalfieldtheory, JournalofHighEnergyPhysics2022(5),152(2022),
doi:10.1007/JHEP05(2022)152.
[50] R. T. Ponnaganti, M. Mambrini and D. Poilblanc, Real-time dynamics of a crit-
ical resonating valence bond spin liquid, Physical Review B 106(19) (2022),
doi:10.1103/physrevb.106.195132.
[51] G. Giudice, G. Giudici, M. Sonner, J. Thoenniss, A. Lerose, D. A. Abanin and L. Piroli,
Temporal entanglement, quasiparticles, and the role of interactions, Phys. Rev. Lett. 128,
220401(2022), doi:10.1103/PhysRevLett.128.220401.
[52] A. Lerose, M. Sonner and D. A. Abanin, Scaling of temporal entangle-
ment in proximity to integrability, Physical Review B 104(3), 035137 (2021),
doi:10.1103/PhysRevB.104.035137.
[53] A.Lerose,M.SonnerandD.A.Abanin, Overcomingtheentanglementbarrierinquantum
many-bodydynamicsviaspace-timeduality, PhysicalReviewB107(6),L060305(2023),
doi:10.1103/PhysRevB.107.L060305.
29

SciPost Physics Codebases Submission
[54] A. Foligno, T. Zhou and B. Bertini, Temporal Entanglement in Chaotic Quantum Circuits,
Physical Review X 13(4), 041008 (2023), doi:10.1103/PhysRevX.13.041008, 2302.
08502.
[55] J. Yao and P. W. Claeys, Temporal Entanglement Profiles in Dual-Unitary Clif-
ford Circuits with Measurements, doi:10.48550/arXiv.2404.14374, Available in:
https://arxiv.org/abs/2404.14374(2024),2404.14374.
[56] T. Prosen and M. Zˇnidariˇc, Is the efficiency of classical simulations of quan-
tum dynamics related to integrability?, Physical Review E 75(1), 015202 (2007),
doi:10.1103/PhysRevE.75.015202.
[57] C. Jonay, D. A. Huse and A. Nahum, Coarse-grained dynamics of operator and state en-
tanglement(2018), 1803.00089.
[58] B. Bertini, P. Kos and T. Prosen, Operator Entanglement in Local Quantum
Circuits I: Chaotic Dual-Unitary Circuits, SciPost Physics 8(4), 067 (2020),
doi:10.21468/SciPostPhys.8.4.067.
[59] M. Heyl, A. Polkovnikov and S. Kehrein, Dynamical Quantum Phase Transitions in
the Transverse Field Ising Model, Physical Review Letters 110(13), 135704 (2013),
doi:10.1103/PhysRevLett.110.135704, 1206.2505.
[60] M.Heyl, Dynamicalquantumphasetransitions: areview, ReportsonProgressinPhysics
81(5),054001(2018), doi:10.1088/1361-6633/aaaf9a.
[61] A. J. Parzygnat, T. Takayanagi, Y. Takiand Z. Wei, Svdentanglement entropy, Journalof
HighEnergyPhysics2023(12)(2023), doi:10.1007/jhep12(2023)123.
[62] P. Caputa, S. Purkayastha, A. Saha and P. Sułkowski, Musings on svd and
pseudo entanglement entropies, Journal of High Energy Physics 2024(11) (2024),
doi:10.1007/jhep11(2024)103.
[63] G. Park, J. Gray and G. K.-L. Chan, Simulating quantum dynamics in two-dimensional
latticeswithtensornetworkinfluencefunctionalbeliefpropagation(2025),2504.07344.
30