# ---- Alcaraz (ANNNI-type) --------------------------------------------------
# H = -Σ_i [ Z_i Z_{i+1} + p λ X_i X_{i+1} + p Z_i Z_{i+2} + λ X_i ]
# Self-dual; p=0 is the integrable TFIM. Field on σx, coupling on σz.

abstract type AbstractAlcarazRecipe <: ExpHRecipe end
struct AlcarazWI <: AbstractAlcarazRecipe end
struct AlcarazWII <: AbstractAlcarazRecipe end
struct AlcarazVD2 <: AbstractAlcarazRecipe end

_alg_string(::AlcarazWI) = "WI"
_alg_string(::AlcarazWII) = "WII"
_alg_string(::AlcarazVD2) = "VD2"

Base.@kwdef mutable struct AlcarazParams <: ModelParams
    lambda::Float64 = 1.0
    p::Float64 = 0.0
    phys_site::Index{Int64} = Index(2, "S=1/2")
end

AlcarazParams(lambda::Number, p::Number) = AlcarazParams(; lambda=Float64(lambda), p=Float64(p))
AlcarazParams(x::AlcarazParams; lambda=x.lambda, p=x.p) = AlcarazParams(; lambda, p, phys_site=x.phys_site)

"""Builds the Alcaraz (ANNNI-type) Hamiltonian as an OpSum."""
function alcaraz_opsum(N::Int, lambda::Number, p::Number)
    os = OpSum()
    for j in 1:(N-1)                 # nearest-neighbour
        os += -1.0, "Z", j, "Z", j + 1
        os += -p * lambda, "X", j, "X", j + 1
    end
    for j in 1:(N-2)                 # next-nearest-neighbour
        os += -p, "Z", j, "Z", j + 2
    end
    for j in 1:N                       # transverse field
        os += -lambda, "X", j
    end
    return os
end



"""Direct U(dt)=exp(-i H dt) MPO for the Alcaraz model. alg = {WI,WII,VD2}"""
function expH_alcaraz(sites::Vector{<:Index}, lambda::Number, p::Number; dt::Number, mpo_alg::String="VD2")
    os = alcaraz_opsum(length(sites), lambda, p)
    return expmpo(os, sites, -im * dt; alg=Algorithm(mpo_alg))
end

# ITransverse: returns the spatial U(dt) MPO that is later rotated into the tMPO
function ITransverse.expH(sites::Vector{<:Index}, mp::AlcarazParams, recipe::AbstractAlcarazRecipe; dt::Number)
    os = alcaraz_opsum(length(sites), mp.lambda, mp.p)
    return expmpo(os, sites, -im * dt; alg=Algorithm(_alg_string(recipe)))
end

# ---- Anisotropic XY --------------------------------------------------------
# H = -Σ_i [ (1+γ)/2 X_i X_{i+1} + (1-γ)/2 Y_i Y_{i+1} + λ Z_i ]
# ITransverse Ising convention (coupling on X, field on Z); γ=1 is the TFIM. Critical at λ=1 for
# any γ>0, Ising class, and the velocity is exactly v = 2γ — which is why we use it to check the
# temporal pipeline really does read a different v at fixed c.

abstract type AbstractXYRecipe <: ExpHRecipe end
struct XYWI <: AbstractXYRecipe end
struct XYWII <: AbstractXYRecipe end
struct XYVD2 <: AbstractXYRecipe end

_alg_string(::XYWI) = "WI"
_alg_string(::XYWII) = "WII"
_alg_string(::XYVD2) = "VD2"

Base.@kwdef mutable struct XYParams <: ModelParams
    lambda::Float64 = 1.0
    gamma::Float64 = 1.0
    phys_site::Index{Int64} = Index(2, "S=1/2")
end

XYParams(lambda::Number, gamma::Number) = XYParams(; lambda=Float64(lambda), gamma=Float64(gamma))
XYParams(x::XYParams; lambda=x.lambda, gamma=x.gamma) = XYParams(; lambda, gamma, phys_site=x.phys_site)

"""Builds the anisotropic XY Hamiltonian as an OpSum (Ising convention: coupling X, field Z)."""
function xy_opsum(N::Int, lambda::Number, gamma::Number)
    os = OpSum()
    for j in 1:(N-1)                 # anisotropic nearest-neighbour coupling
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