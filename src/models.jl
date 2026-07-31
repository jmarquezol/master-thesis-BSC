# ── Alcaraz (ANNNI-type) ────────────────────────────────────────────────────────────────────
# H = -Σ_i [ Z_i Z_{i+1} + p λ X_i X_{i+1} + p Z_i Z_{i+2} + λ X_i ]
# Self-dual; p=0 is the integrable TFIM. Field on σx, coupling on σz.

abstract type AbstractAlcarazRecipe <: ExpHRecipe end
struct AlcarazWI  <: AbstractAlcarazRecipe end
struct AlcarazWII <: AbstractAlcarazRecipe end
struct AlcarazVD2 <: AbstractAlcarazRecipe end

_alg_string(::AlcarazWI)  = "WI"
_alg_string(::AlcarazWII) = "WII"
_alg_string(::AlcarazVD2) = "VD2"

Base.@kwdef mutable struct AlcarazParams <: ModelParams
    lambda::Float64 = 1.0
    p::Float64      = 0.0
    phys_site::Index{Int64} = Index(2, "S=1/2")
end

AlcarazParams(lambda::Number, p::Number) = AlcarazParams(; lambda=Float64(lambda), p=Float64(p))
AlcarazParams(x::AlcarazParams; lambda=x.lambda, p=x.p) = AlcarazParams(; lambda, p, phys_site=x.phys_site)

"""Builds the Alcaraz (ANNNI-type) Hamiltonian as an OpSum."""
function alcaraz_opsum(N::Int, lambda::Number, p::Number)
    os = OpSum()
    for j in 1:(N - 1)                 # nearest-neighbour
        os += -1.0,        "Z", j, "Z", j + 1
        os += -p * lambda, "X", j, "X", j + 1
    end
    for j in 1:(N - 2)                 # next-nearest-neighbour
        os += -p, "Z", j, "Z", j + 2
    end
    for j in 1:N                       # transverse field
        os += -lambda, "X", j
    end
    return os
end

 

"""Direct U(dt)=exp(-i H dt) MPO for the Alcaraz model. alg = {WI,WII,VD2}"""
function expH_alcaraz(sites::Vector{<:Index}, lambda::Number, p::Number; dt::Number, mpo_alg::String="VD2")
    os  = alcaraz_opsum(length(sites), lambda, p)
    return expmpo(os, sites, -im * dt; alg=Algorithm(mpo_alg))
end

# ITransverse: returns the spatial U(dt) MPO that is later rotated into the tMPO
function ITransverse.expH(sites::Vector{<:Index}, mp::AlcarazParams, recipe::AbstractAlcarazRecipe; dt::Number)
    os = alcaraz_opsum(length(sites), mp.lambda, mp.p)
    return expmpo(os, sites, -im * dt; alg=Algorithm(_alg_string(recipe)))
end

# ── Anisotropic XY ──────────────────────────────────────────────────────────────────────────
# H = -Σ_i [ (1+γ)/2 X_i X_{i+1} + (1-γ)/2 Y_i Y_{i+1} + λ Z_i ]
# ITransverse Ising convention (coupling on X, field on Z); γ=1 is the TFIM. Critical at λ=1 for
# any γ>0, Ising class, and the velocity is exactly v = 2γ — which is why we use it to check the
# temporal pipeline really does read a different v at fixed c.

abstract type AbstractXYRecipe <: ExpHRecipe end
struct XYWI  <: AbstractXYRecipe end
struct XYWII <: AbstractXYRecipe end
struct XYVD2 <: AbstractXYRecipe end

_alg_string(::XYWI)  = "WI"
_alg_string(::XYWII) = "WII"
_alg_string(::XYVD2) = "VD2"

Base.@kwdef mutable struct XYParams <: ModelParams
    lambda::Float64 = 1.0
    gamma::Float64  = 1.0
    phys_site::Index{Int64} = Index(2, "S=1/2")
end

XYParams(lambda::Number, gamma::Number) = XYParams(; lambda=Float64(lambda), gamma=Float64(gamma))
XYParams(x::XYParams; lambda=x.lambda, gamma=x.gamma) = XYParams(; lambda, gamma, phys_site=x.phys_site)

"""Builds the anisotropic XY Hamiltonian as an OpSum (Ising convention: coupling X, field Z)."""
function xy_opsum(N::Int, lambda::Number, gamma::Number)
    os = OpSum()
    for j in 1:(N - 1)                 # anisotropic nearest-neighbour coupling
        os += -(1 + gamma) / 2, "X", j, "X", j + 1
        os += -(1 - gamma) / 2, "Y", j, "Y", j + 1
    end
    for j in 1:N                       # transverse field
        os += -lambda, "Z", j
    end
    return os
end

"""Direct U(dt)=exp(-i H dt) MPO for the XY model. alg = {WI,WII,VD2}"""
function expH_xy(sites::Vector{<:Index}, lambda::Number, gamma::Number; dt::Number, mpo_alg::String="VD2")
    os = xy_opsum(length(sites), lambda, gamma)
    return expmpo(os, sites, -im * dt; alg=Algorithm(mpo_alg))
end

function ITransverse.expH(sites::Vector{<:Index}, mp::XYParams, recipe::AbstractXYRecipe; dt::Number)
    os = xy_opsum(length(sites), mp.lambda, mp.gamma)
    return expmpo(os, sites, -im * dt; alg=Algorithm(_alg_string(recipe)))
end

# ── Tricritical Ising (optional variant) ────────────────────────────────────────────────────
# H = -Σ X_i - Σ Z_i Z_{i+1} + λ Σ ( Z_i Z_{i+1} X_{i+2} + X_i Z_{i+1} Z_{i+2} )

abstract type AbstractTricriticalRecipe <: ExpHRecipe end
struct TricriticalWI  <: AbstractTricriticalRecipe end
struct TricriticalWII <: AbstractTricriticalRecipe end
struct TricriticalVD2 <: AbstractTricriticalRecipe end

_alg_string(::TricriticalWI)  = "WI"
_alg_string(::TricriticalWII) = "WII"
_alg_string(::TricriticalVD2) = "VD2"

Base.@kwdef mutable struct TricriticalParams <: ModelParams
    lambda::Float64 = 0.0
    phys_site::Index{Int64} = Index(2, "S=1/2")
end

TricriticalParams(lambda::Number) = TricriticalParams(; lambda=Float64(lambda))
TricriticalParams(x::TricriticalParams; lambda=x.lambda) = TricriticalParams(; lambda, phys_site=x.phys_site)

"""Builds the tricritical-Ising Hamiltonian (with 3-body terms) as an OpSum."""
function tricritical_opsum(N::Int, lambda::Number)
    os = OpSum()
    for j in 1:N                       # transverse field
        os += -1.0, "X", j
    end
    for j in 1:(N - 1)                 # nearest-neighbour
        os += -1.0, "Z", j, "Z", j + 1
    end
    for j in 1:(N - 2)                 # 3-body terms
        os += lambda, "Z", j, "Z", j + 1, "X", j + 2
        os += lambda, "X", j, "Z", j + 1, "Z", j + 2
    end
    return os
end

"""Direct U(dt) MPO for the tricritical model (Schrödinger pipeline)."""
function expH_tricritical(sites::Vector{<:Index}, lambda::Number; dt::Number, mpo_alg::String="VD2")
    os  = tricritical_opsum(length(sites), lambda)
    return expmpo(os, sites, -im * dt; alg=Algorithm(mpo_alg))
end

function ITransverse.expH(sites::Vector{<:Index}, mp::TricriticalParams, recipe::AbstractTricriticalRecipe; dt::Number)
    os = tricritical_opsum(length(sites), mp.lambda)
    return expmpo(os, sites, -im * dt; alg=Algorithm(_alg_string(recipe)))
end

# ── XXZ, asymmetric exp-MPO (VD2) path ──────────────────────────────────────────────────────
# Reuses ITransverse's XXZParams (J_XY, J_ZZ, hz): H = -( J_XY (XX+YY + J_ZZ ZZ) + 2 hz Z ), so
# H_Δ = Σ_x [½(S+S- + S-S+) + Δ Sz Sz] is XXZParams(-1.0, Δ, 0.0). The overall sign of H does not
# matter for |L| or the temporal entropy. ITransverse ships a symmetric SymSVD builder; this is
# the 2nd-order asymmetric one, for the comparison against Alcaraz.

abstract type AbstractXXZRecipe <: ExpHRecipe end
struct XXZWI  <: AbstractXXZRecipe end
struct XXZWII <: AbstractXXZRecipe end
struct XXZVD2 <: AbstractXXZRecipe end

_alg_string(::XXZWI)  = "WI"
_alg_string(::XXZWII) = "WII"
_alg_string(::XXZVD2) = "VD2"

"""Builds the XXZ Hamiltonian as an OpSum (ITransverse convention H=-(J_XY(XX+YY+Δ ZZ)+2hZ))."""
function xxz_opsum(N::Int, J_XY::Number, J_ZZ::Number, hz::Number=0.0)
    os = OpSum()
    for j in 1:(N - 1)
        os += -J_XY / 2,    "S+", j, "S-", j + 1
        os += -J_XY / 2,    "S-", j, "S+", j + 1
        if abs(J_ZZ) > 1e-12
            os += -J_XY * J_ZZ, "Sz", j, "Sz", j + 1
        end
    end
    if abs(hz) > 1e-12
        for j in 1:N
            os += -2 * hz, "Sz", j
        end
    end
    return os
end

"""Direct U(dt)=exp(-i H dt) MPO for the XXZ model (Schrödinger pipeline). alg = {WI,WII,VD2}"""
function expH_xxz(sites::Vector{<:Index}, J_XY::Number, J_ZZ::Number; dt::Number, mpo_alg::String="VD2", hz::Number=0.0)
    os = xxz_opsum(length(sites), J_XY, J_ZZ, hz)
    return expmpo(os, sites, -im * dt; alg=Algorithm(mpo_alg))
end

function ITransverse.expH(sites::Vector{<:Index}, mp::XXZParams, recipe::AbstractXXZRecipe; dt::Number)
    os = xxz_opsum(length(sites), mp.J_XY, mp.J_ZZ, mp.hz)
    return expmpo(os, sites, -im * dt; alg=Algorithm(_alg_string(recipe)))
end

# ── XXZ Néel quench, in the sublattice-rotated frame ────────────────────────────────────────
# We want to quench |↑↓↑↓…⟩ with H_Δ = Σ_x [½(S+S- + S-S+) + Δ Sz Sz]. The rotation
# R = ∏_{x even} exp(iπ Sx) sends |Néel⟩ → |↑↑↑…⟩ and H_Δ → H'_Δ, with S+S- → S+S+ and −Δ on the
# ZZ term. R is a product of single-site unitaries, so the echo and the whole temporal-entropy
# structure are unchanged — we evolve the uniform |↑⟩ under H'_Δ and keep the single-site
# machinery. Note `Delta` stores the PHYSICAL +Δ; the sign flip in the rotated frame is forced,
# not a typo (a Néel→ferromagnet rotation is a π-rotation in the XY plane, which flips Sz).
# Checked against a direct Néel TDVP echo to 4 digits. Critical for |Δ| ≤ 1.

abstract type AbstractXXZNeelRecipe <: ExpHRecipe end
struct XXZNeelWI  <: AbstractXXZNeelRecipe end
struct XXZNeelWII <: AbstractXXZNeelRecipe end
struct XXZNeelVD2 <: AbstractXXZNeelRecipe end

_alg_string(::XXZNeelWI)  = "WI"
_alg_string(::XXZNeelWII) = "WII"
_alg_string(::XXZNeelVD2) = "VD2"

Base.@kwdef mutable struct XXZNeelParams <: ModelParams
    Delta::Float64 = 0.5                      # PHYSICAL XXZ anisotropy (supervisor's +Δ SzSz)
    phys_site::Index{Int64} = Index(2, "S=1/2")
end

XXZNeelParams(Δ::Number) = XXZNeelParams(; Delta=Float64(Δ))
XXZNeelParams(x::XXZNeelParams; Delta=x.Delta) = XXZNeelParams(; Delta, phys_site=x.phys_site)

"""Sublattice-rotated XXZ Hamiltonian (Néel→uniform |↑⟩ frame): Σ ½(S+S+ + S-S-) − Δ SzSz.
   `Delta` is the physical (+Δ) anisotropy; the −Δ here is the rotation-induced ZZ sign flip."""
function xxz_neel_opsum(N::Int, Delta::Number)
    os = OpSum()
    for j in 1:(N - 1)
        os += 0.5,      "S+", j, "S+", j + 1
        os += 0.5,      "S-", j, "S-", j + 1
        os += -Delta,   "Sz", j, "Sz", j + 1
    end
    return os
end

"""Direct U(dt) MPO for the rotated XXZ-Néel model (Schrödinger pipeline)."""
function expH_xxz_neel(sites::Vector{<:Index}, Delta::Number; dt::Number, mpo_alg::String="VD2")
    os = xxz_neel_opsum(length(sites), Delta)
    return expmpo(os, sites, -im * dt; alg=Algorithm(mpo_alg))
end

function ITransverse.expH(sites::Vector{<:Index}, mp::XXZNeelParams, recipe::AbstractXXZNeelRecipe; dt::Number)
    os = xxz_neel_opsum(length(sites), mp.Delta)
    return expmpo(os, sites, -im * dt; alg=Algorithm(_alg_string(recipe)))
end

# ── symmetric (Murg-type) propagator for the rotated XXZ-Néel chain ─────────────────────────
# In Pauli form H'_Δ = Σ_j ¼(σx σx − σy σy − Δ σz σz). Each layer Σ_j σᵃ_j σᵃ_{j+1} commutes with
# itself, so its exponential is an exact bond-2 left-right-symmetric MPO (the Murg cos/sin split,
# as in ITransverse's expXX_murg but for any Pauli a). That is what opens the symmetric route
# (powermethod_sym + Takagi) for XXZ, which the asymmetry experiment needed.
#
# Two kernels, both reflection-symmetric — the palindrome buys Trotter order, not symmetry:
#   order=2 (default): e^{ZZ/2}e^{YY/2}e^{XX}e^{YY/2}e^{ZZ/2}, 2nd order, temporal dim d_t=32.
#   order=1: e^{ZZ}e^{YY}e^{XX}, 1st order but d_t=8, so ~16x cheaper per apply. Validate the
#     O(dt) error at your dt before using it (NB12 §1).

struct XXZNeelMurg <: AbstractXXZNeelRecipe
    order::Int
end
XXZNeelMurg() = XXZNeelMurg(2)

"""Exact symmetric bond-2 MPO of exp(+i*Jdt*Σ_j σᵃ_j σᵃ_{j+1}), a ∈ {"X","Y","Z"} (Pauli ops).
   Transcription of ITransverse's `expXX_murg` cos/sin splitting, generalized to any Pauli."""
function exp2site_murg(sites::Vector{<:Index}, Jdt::Number, opname::String)
    N = length(sites)
    U = MPO(N)
    links = [Index(2, "Link,l=$(n-1)") for n in 1:(N + 1)]
    for n in 1:N
        ll, rl = dag(links[n]), links[n + 1]
        I = op(sites, "Id", n)
        A = op(sites, opname, n)
        if n == 1
            U[n]  = onehot(rl => 1) * sqrt(cos(Jdt)) * I
            U[n] += onehot(rl => 2) * sqrt(im * sin(Jdt)) * A
        elseif n == N
            U[n]  = onehot(ll => 1) * sqrt(cos(Jdt)) * I
            U[n] += onehot(ll => 2) * sqrt(im * sin(Jdt)) * A
        else
            U[n]  = onehot(ll => 1, rl => 1) * cos(Jdt) * I
            U[n] += onehot(ll => 1, rl => 2) * sqrt(im * sin(Jdt)) * sqrt(cos(Jdt)) * A
            U[n] += onehot(ll => 2, rl => 1) * sqrt(im * sin(Jdt)) * sqrt(cos(Jdt)) * A
            U[n] += onehot(ll => 2, rl => 2) * im * sin(Jdt) * I
        end
    end
    return U
end

"""U(dt) = exp(−i dt H'_Δ) via a Murg sandwich of the mutually-commuting XX/YY/ZZ layers.
   `order=2` (default): the palindrome e^{ZZ/2}e^{YY/2}e^{XX}e^{YY/2}e^{ZZ/2}, 2nd order in dt,
   d_t=32. `order=1`: the single sandwich e^{ZZ}e^{YY}e^{XX} (full steps), 1st order in dt, d_t=8.
   Layer coefficients (exp(−i dt J Σσσ) = exp(+i(−J dt)Σσσ), spin→Pauli factor ¼):
   XX: J=+¼ → Jdt = −dt/4 (full) or −dt/4 (order=1, same — XX is never split);
   YY: J=−¼ → Jdt = +dt/8 (half, order=2) or +dt/4 (full, order=1);
   ZZ: J=−Δ/4 → Jdt = +Δdt/8 (half, order=2) or +Δdt/4 (full, order=1)."""
function expH_xxz_neel_murg(sites::Vector{<:Index}, Delta::Number; dt::Number, order::Int=2)
    Uxx = exp2site_murg(sites, -dt / 4, "X")
    if order == 1
        Uzz = exp2site_murg(sites, Delta * dt / 4, "Z")
        Uyy = exp2site_murg(sites, dt / 4,         "Y")
        U = applyn(Uyy, Uzz)            # build e^{ZZ}·e^{YY}·e^{XX}
        U = applyn(Uxx, U)
        return U
    elseif order == 2
        Uzz2 = exp2site_murg(sites, Delta * dt / 8, "Z")
        Uyy2 = exp2site_murg(sites, dt / 8,         "Y")
        U = applyn(Uyy2, Uzz2)          # build the palindrome Uzz2·Uyy2·Uxx·Uyy2·Uzz2
        U = applyn(Uxx,  U)
        U = applyn(Uyy2, U)
        U = applyn(Uzz2, U)
        return U
    else
        error("expH_xxz_neel_murg: order must be 1 or 2, got $order")
    end
end

ITransverse.expH(sites::Vector{<:Index}, mp::XXZNeelParams, recipe::XXZNeelMurg; dt::Number) =
    expH_xxz_neel_murg(sites, mp.Delta; dt=dt, order=recipe.order)
