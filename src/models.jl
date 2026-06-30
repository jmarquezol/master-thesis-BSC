# ══════════════════════════════════════════════════════════════════════════════
# ALCARAZ (ANNNI-type) MODEL
#   H = -Σ_i [ Z_i Z_{i+1} + p λ X_i X_{i+1} + p Z_i Z_{i+2} + λ X_i ]
#   self-dual; p=0 recovers the integrable TFIM. Field on sigma_x, coupling on sigma_z
# ══════════════════════════════════════════════════════════════════════════════

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

# ══════════════════════════════════════════════════════════════════════════════
# TRICRITICAL-ISING MODEL (optional variant)
#   H = -Σ X_i - Σ Z_i Z_{i+1} + λ Σ ( Z_i Z_{i+1} X_{i+2} + X_i Z_{i+1} Z_{i+2} )
# ══════════════════════════════════════════════════════════════════════════════

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

# ══════════════════════════════════════════════════════════════════════════════
# XXZ MODEL — our asymmetric exp-MPO (VD2) path
#   Reuses ITransverse's XXZParams (J_XY, J_ZZ, hz) with convention
#     H = -( J_XY (XX+YY + J_ZZ ZZ) + 2 hz Z ),   J_ZZ = Δ.
#   The supervisor's  H_Δ = Σ_x [½(S+S- + S-S+) + Δ Sz Sz]  ⇒  XXZParams(-1.0, Δ, 0.0)
#   (overall sign of H is immaterial to |L| and the temporal entropy).
#   ITransverse already provides the SYMMETRIC SymSVD builder; THIS adds the 2nd-order
#   asymmetric VD2 builder (like Alcaraz) for the gap-vs-Alcaraz comparison and a
#   higher-Trotter-order cross-check.
# ══════════════════════════════════════════════════════════════════════════════

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

# ══════════════════════════════════════════════════════════════════════════════
# XXZ NÉEL QUENCH (sublattice-rotated frame) — the model we actually time-evolve
#   Physical setup: quench the NÉEL state |↑↓↑↓…⟩ with the critical XXZ Hamiltonian
#     H_Δ = Σ_x [ ½(S+_x S-_{x+1} + S-_x S+_{x+1}) + Δ Sz_x Sz_{x+1} ]   (supervisor's Eq. 2, +Δ)
#   The single-site sublattice rotation R = ∏_{x even} exp(iπ Sx_x) maps
#     |Néel⟩ → |↑↑↑…⟩ (uniform)  and
#     H_Δ → H'_Δ = Σ_x [ ½(S+_x S+_{x+1} + S-_x S-_{x+1}) − Δ Sz_x Sz_{x+1} ].
#   R is a product of single-site unitaries ⇒ the Loschmidt echo and the entire temporal-
#   entropy structure are IDENTICAL, so we evolve the UNIFORM |↑⟩ state under H'_Δ and reuse
#   the single-site transverse machinery (no 2-site unit cell needed).
#   IMPORTANT: `Delta` is the PHYSICAL XXZ anisotropy (supervisor's +Δ convention); the ZZ
#   sign necessarily flips to −Δ in the rotated frame (any Néel→ferromagnet rotation is a
#   π-rotation in the XY plane, which flips Sz). Verified: |↑⟩-under-H'_Δ echo == direct
#   Néel-under-H_Δ TDVP echo to 4 digits. Critical regime |Δ| ≤ 1 (c = 1 Luttinger liquid).
# ══════════════════════════════════════════════════════════════════════════════

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
