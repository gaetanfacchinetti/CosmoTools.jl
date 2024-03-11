##################################################################################
# This file is part of CosmoTools.jl
#
# Copyright (c) 2024, Gaétan Facchinetti
#
# CosmoTools.jl is free software: you can redistribute it and/or modify it 
# under the terms of the GNU General Public License as published by 
# the Free Software Foundation, either version 3 of the License, or any 
# later version. CosmoTools.jl is distributed in the hope that it will be useful, 
# but WITHOUT ANY WARRANTY; without even the implied warranty of 
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. 
# See the GNU General Public License for more details.
#
# You should have received a copy of the GNU 
# General Public License along with 21cmCAST. 
# If not, see <https://www.gnu.org/licenses/>.
##################################################################################


export Halo, nfwProfile, αβγProfile, HaloProfile, coreProfile, plummerProfile
export halo_from_mΔ_and_cΔ
export mΔ_from_ρs_and_rs, mΔ, rΔ_from_ρs_and_rs, rΔ, cΔ_from_ρs, cΔ, ρ_halo, μ_halo, m_halo
export velocity_dispersion, gravitational_potential, escape_velocity, orbital_frequency
export gravitational_potential_kms

abstract type HaloProfile{T<:Real} end

################################################
# Dimensionless Halo

struct αβγProfile{T<:Real} <: HaloProfile{T}
    name::String
    α::T
    β::T
    γ::T
end

# definition of length and iterator on our struct
# allows to use f.(x, y) where y is of type HaloProfile
Base.length(::HaloProfile) = 1
Base.iterate(iter::HaloProfile) = (iter, nothing)
Base.iterate(::HaloProfile, state::Nothing) = nothing

# overrinding print function
Base.show(io::IO, hp::αβγProfile{<:Real}) = print(io, "αβγProfile: α = " * string(hp.α) * ", β = " * string(hp.β) * ", γ = " * string(hp.γ))


# constructor method of αβγProfile
αβγProfile(name::String, α::Real, β::Real, γ::Real) = αβγProfile(name, promote(α, β, γ)...)

## definition of densities and mass profiles
const nfwProfile::αβγProfile = αβγProfile("NFW_Profile", 1, 3, 1)
const coreProfile::αβγProfile = αβγProfile("Core_Profile", 1, 3, 0)
const plummerProfile::αβγProfile = αβγProfile("Plummer_Profile", 2, 5, 0)

# Density profile and mass
ρ_halo(x::Real, p::αβγProfile = nfwProfile) = x^(-p.γ) * (1+x^p.α)^(-(p.β - p.γ)/p.α)
μ_halo(x::Real, p::αβγProfile = nfwProfile) = HypergeometricFunctions._₂F₁((3 - p.γ)/p.α, (p.β - p.γ)/p.α, (3 + p.α - p.γ)/p.α, -x^p.α) * x^(3-p.γ) / (3-p.γ)

function gravitational_potential(x::Real, xt::Real, p::αβγProfile = nfwProfile) 
    (p == nfwProfile) && return (log(1+x)/x -log(1+xt)/xt)
    return - quadgk(xp -> μ_halo(xp, p) / xp^2, x, xt, rtol=1e-3)[1] 
end

# s- and p-wave luminosity of the halo
λ0_halo(x::Real, p::αβγProfile = nfwProfile) = HypergeometricFunctions._₂F₁((3 - 2*p.γ)/p.α, 2*(p.β - p.γ)/p.α, (3 + p.α - 2*p.γ)/p.α, -x^p.α) * x^(3-2*p.γ) / (3-2*p.γ)
λ1_halo(x::Real, p::αβγProfile = nfwProfile) = 1.0 # TO DO

################################################
#
# Dimensionfull Halo

# Whenever possible do all the computation using the dimensionless profile
# Otherwise unit conversion slightly slow down the code

## Definition of relationships between concentration, mass, scale density and scale radius
function cΔ_from_ρs(ρs::Real, hp::HaloProfile{<:Real}, Δ::Real, ρ_ref::Real)
    g(c::Real) = c^3 / μ_halo(c, hp) - 3 * ρs / Δ / ρ_ref
    Roots.find_zero(g, (1e-10, 1e+10), Roots.Bisection()) 
end

cΔ_from_ρs(ρs::Real, hp::HaloProfile{<:Real} = nfwProfile, Δ::Real=200, cosmo::Cosmology = plnack18) = cΔ_from_ρs(ρs, hp, Δ, cosmo.bkg.ρ_c0)

mΔ_from_ρs_and_rs(ρs::Real, rs::Real, hp::HaloProfile{<:Real} = nfwProfile, Δ::Real = 200, ρ_ref::Real = planck18_bkg.ρ_c0) = 4 * pi * ρs * rs^3 * μ_halo(cΔ_from_ρs(ρs, hp, Δ, ρ_ref), hp)
ρs_from_cΔ(cΔ::Real, hp::HaloProfile{<:Real} = nfwProfile , Δ::Real = 200, ρ_ref::Real = planck18_bkg.ρ_c0) = Δ * ρ_ref  / 3 * cΔ^3 / μ_halo(cΔ, hp) 
rs_from_cΔ_and_mΔ(cΔ::Real, mΔ::Real, Δ::Real = 200, ρ_ref::Real = planck18_bkg.ρ_c0) =  (3 * mΔ / (4 * π * Δ * ρ_ref))^(1 // 3) / cΔ 
rΔ_from_ρs_and_rs(ρs::Real, rs::Real, hp::HaloProfile{<:Real} = nfwProfile, Δ::Real = 200, ρ_ref::Real = planck18_bkg.ρ_c0) = (3 * mΔ_from_ρs_and_rs(ρs, rs, hp, Δ, ρ_ref) / (4*π*Δ*ρ_ref))^(1//3)

struct Halo{T<:Real}
    hp::HaloProfile
    ρs::T
    rs::T
end

Base.length(::Halo) = 1
Base.iterate(iter::Halo) = (iter, nothing)
Base.iterate(::Halo, state::Nothing) = nothing

Base.show(io::IO, h::Halo{<:Real}) = print(io, "Halo: \n  - " * string(h.hp) * "\n  - ρs = " * string(h.ρs) * " Msun/Mpc^3, rs = " * string(h.rs) * " Mpc \n  - m200 (planck18) = " * string(mΔ(h, 200, planck18)) * " Msun, c200 (planck18) = " * string(cΔ(h, 200, planck18)))

Halo(hp::HaloProfile, ρs::Real, rs::Real) = Halo(hp, promote(ρs, rs)...)

function halo_from_mΔ_and_cΔ(hp::HaloProfile, mΔ::Real, cΔ::Real;  Δ::Real = 200, ρ_ref::Real = planck18_bkg.ρ_c0)
    ρs = ρs_from_cΔ(cΔ, hp, Δ, ρ_ref)
    rs = rs_from_cΔ_and_mΔ(cΔ, mΔ, Δ, ρ_ref)
    return Halo(hp, ρs, rs)
end

ρ_halo(r::Real, h::Halo{<:Real}) = h.ρs * ρ_halo(r/h.rs, h.hp)
m_halo(r::Real, h::Halo{<:Real}) = 4.0 * π * h.ρs * h.rs^3 * μ_halo(r/h.rs, h.hp)
mΔ(h::Halo{<:Real}, Δ::Real, ρ_ref::Real) = mΔ_from_ρs_and_rs(h.ρs, h.rs, h.hp, Δ, ρ_ref)
mΔ(h::Halo{<:Real}, Δ::Real = 200, cosmo::Cosmology = planck18) = mΔ(h, Δ, cosmo.bkg.ρ_c0)
cΔ(h::Halo{<:Real}, Δ::Real, ρ_ref::Real) = cΔ_from_ρs(h.ρs, h.hp, Δ, ρ_ref)
cΔ(h::Halo{<:Real}, Δ::Real = 200, cosmo::Cosmology = planck18) = cΔ(h, Δ, cosmo.bkg.ρ_c0)
rΔ(h::Halo{<:Real}, Δ::Real, ρ_ref::Real) = rΔ_from_ρs_and_rs(h.ρs, h.rs, h.hp, Δ, ρ_ref)
rΔ(h::Halo{<:Real}, Δ::Real = 200, cosmo::Cosmology = planck18) = rΔ(h, Δ, cosmo.bkg.ρ_c0)

ρs(h::Halo) = h.ρs * ρ_0
rs(h::Halo) = h.rs * r_0

################################################
# Velocity of the particles inside the halo

""" 1D Jeans velocity dispersion in (Mpc / s) """
velocity_dispersion(r::Real, rt::Real, h::Halo) = sqrt(G_NEWTON / ρ_halo(r, h) *  quadgk(rp -> ρ_halo(rp, h) * m_halo(rp, h)/rp^2, r, rt, rtol=1e-5)[1])

""" 1D Jeans velocity dispersion in (km / s) """
velocity_dispersion_kms(r::Real, rt::Real, h::Halo) = velocity_dispersion(r, rt, h) * MPC_TO_KM

""" gravitational potential in (Mpc / s)^2 """
gravitational_potential(r::Real, rt::Real, h::Halo) = - 4 * π * h.ρs * h.rs^2 * G_NEWTON * gravitational_potential(r/h.rs, rt/h.rs, h.hp) 

""" gravitational potential in (km / s)^2 """
gravitational_potential_kms(r::Real, rt::Real, h::Halo) = gravitational_potential(r, rt, h) * (MPC_TO_KM^2)

""" escape velocity in (Mpc / s) """
escape_velocity(r::Real, rt::Real, h::Halo) = sqrt(2*abs(gravitational_potential(r, rt, h)))

""" escape velocity in (km / s) """
escape_velocity_kms(r::Real, rt::Real, h::Halo) = sqrt(2*abs(gravitational_potential_kms(r, rt, h))) * MPC_TO_KM

""" circular velocity in (Mpc / s)"""
circular_velocity(r::Real, h::Halo) = sqrt(G_NEWTON * m_halo(r, h) / r)

""" circular velocity in (km / s)"""
circular_velocity_kms(r::Real, h::Halo) = sqrt(G_NEWTON * m_halo(r, h) / r) * MPC_TO_KM

""" orbital frequency in (1 / s) """
orbital_frequency(r::Real, h::Halo) = circular_velocity(r, h) / (2*π*r)



#######################################
# Relation between concentration and mass
# Considering relations derived in the litterature

export MassConcentrationModel, SCP12, MDvdB08, median_concentration, PKCBRP12

abstract type MassConcentrationModel end
abstract type SCP12 <: MassConcentrationModel end
abstract type MDvdB08 <: MassConcentrationModel end
abstract type PKCBRP12 <: MassConcentrationModel end

function median_concentration(m200::Real, z::Real, cosmo::Cosmology, ::Type{SCP12})

    (z != 0) && Throw(ArgumentError("The mass-concentration model SCP12 is not defined for z > 0"))

    cn::Vector = [37.5153, -1.5093, 1.636e-2, 3.66e-4, -2.89237e-5, 5.32e-7]
    m200_min = (m200 > 7.24e-10) ? m200 : 7.24e-10

    return sum(cn .* log(m200_min * cosmo.bkg.h).^(0:5))
end


function median_concentration(m200::Real, z::Real, cosmo::Cosmology, ::Type{MDvdB08})
    
    # parameters of the model
    K200 = 3.6
    F = 0.01

    zc = 10^Roots.find_zero(logz -> σ_mps_M(F * m200, TopHatWindow, cosmology = cosmo) * growth_factor_Carroll(10^(logz), cosmo.bkg) / growth_factor_Carroll(0.0, cosmo.bkg) - 1.686, (-6, 10), Roots.Bisection())
    return K200 * (ρ_critical(zc, cosmo.bkg)/ ρ_critical(z, cosmo.bkg))^(1/3)
end


function median_concentration(m200::Real, z::Real, cosmo::Cosmology, ::Type{PKCBRP12})
    
    x = (cosmo.bkg.Ω_Λ0 / cosmo.bkg.Ω_m0)^(1/3)/(1+z)
    
    func_C(sigp::Real) = 2.881*( (sigp/1.257)^1.022 +1  )*exp(0.060/sigp^2)
    func_cmin(x::Real) = 3.681 + (5.033 - 3.681)*(atan(6.948 * (x - 0.424))/π + 0.5)
    func_sigm1min(x::Real) = 1.047 + (1.646 - 1.047)*(atan(7.386 * (x - 0.526))/π + 0.5)
    
    func_B0(x::Real) = func_cmin(x)/func_cmin(1.393)
    func_B1(x::Real) = func_sigm1min(x) / func_sigm1min(1.393)

    norm_D = growth_factor_Carroll(z, cosmo.bkg) /  growth_factor_Carroll(0, cosmo.bkg) 
    sigp   = func_B1(x) * σ_mps_M(m200, TopHatWindow, cosmology = cosmo) * norm_D


    return func_B0(x) * func_C(sigp)
end


median_concentration(m200::Real, ::Type{T} = SCP12; z::Real = 0, cosmo::Cosmology=planck18) where {T <: MassConcentrationModel} = median_concentration(m200, z, cosmo, T)

#######################################


