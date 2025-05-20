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
export halo_from_mΔ_and_cΔ, halo_from_mΔ
export mΔ_from_ρs_and_rs, mΔ, rΔ_from_ρs_and_rs, rΔ, cΔ_from_ρs, cΔ, ρ_halo, μ_halo, m_halo
export rs_from_cΔ_and_mΔ, ρs_from_cΔ
export velocity_dispersion, gravitational_potential, escape_velocity, orbital_frequency, circular_velocity
export gravitational_potential_kms
export get_halo_profile_type

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

get_halo_profile_type(::HaloProfile{T}) where {T<:Real} = T

# overrinding print function
Base.show(io::IO, hp::αβγProfile{<:Real}) = print(io, "αβγProfile: α = " * string(hp.α) * ", β = " * string(hp.β) * ", γ = " * string(hp.γ))

# constructor method of αβγProfile
αβγProfile(name::String, α::Real, β::Real, γ::Real) = αβγProfile(name, promote(α, β, γ)...)

## definition of densities and mass profiles
const nfwProfile::αβγProfile{Int} = αβγProfile("NFW_Profile", 1, 3, 1)
const coreProfile::αβγProfile{Int} = αβγProfile("Core_Profile", 1, 3, 0)
const plummerProfile::αβγProfile{Int} = αβγProfile("Plummer_Profile", 2, 5, 0)

# Density profile and mass
ρ_halo(x::Real, p::αβγProfile = nfwProfile) = x^(-p.γ) * (1+x^p.α)^(-(p.β - p.γ)/p.α)

function μ_halo(x::Real, p::αβγProfile = nfwProfile) 
    (p == nfwProfile) && return (x < 1e-3 ? (x^2/2 - (2 * x^3)/3 + (3 * x^4)/4 - (4 * x^5)/5) : (log(1+x) - x/(1+x)))
    return HypergeometricFunctions._₂F₁((3 - p.γ)/p.α, (p.β - p.γ)/p.α, (3 + p.α - p.γ)/p.α, -x^p.α) * x^(3-p.γ) / (3-p.γ)
end

# dimensionless properties of the halo
# use analytical expression when possible
function ψ_halo(x::Real, xt::Real, p::αβγProfile = nfwProfile) 
    (p == nfwProfile) && ( (x > 0) && return (log(1+x)/x -log(1+xt)/xt))  && ( (x == 0) && return (1 -log(1+xt)/xt))
    return - QuadGK.quadgk(lnxp -> μ_halo(exp(lnxp), p) / exp(lnxp), x, xt, rtol=1e-5)[1] 
end


#-------------
#
# Usefull functions for Eddington's inversion
# 

# first derivative of the halo mass density
der_1_ρ_halo(x::Real, p::αβγProfile = nfwProfile) = -x^(-1-p.γ) * (1+x^p.α)^((p.α + p.β - p.γ)/p.α) * (p.β * x^(p.α) + p.γ)

# second derivative of the halo mass density
der_2_ρ_halo(x::Real, p::αβγProfile = nfwProfile) = x^(-2-p.γ) * (1+x^p.α)^((-2 * p.α-p.β+p.γ)/p.α) * (x^(2 * p.α) * p.β *  (1+p.β)+p.γ+p.γ^2+x^p.α * (p.β-p.α * p.β+p.γ+p.α * p.γ+ 2 * p.β * p.γ))

# derivative of the gravitational potential with respect to the mass
der_1_ψ_halo(x::Real, p::HaloProfile = nfwProfile) = -μ_halo(x, p) / x^2

# second derivative of the density profile with respect to the gravitational potential
function der_2_ρ_halo_vs_ψ(x::Real, p::HaloProfile = nfwProfile)
  m = μ_halo(x, p)
  return (2 * x^3) / m^2 * ( (1 - x^3 * ρ_halo(x, p)  / (2 * m)) * der_1_ρ_halo(x, p) + (x / 2) * der_2_ρ_halo(x, p)  );
end 

# for the Eddington inversion method we define φ = ψ/E
# for a given value of E and φ we must get the associated position x
function x_vs_φ(φ::Real, E::Real, x0::Real, p::HaloProfile = nfwProfile)
  _to_solve(x::Real) = ψ_halo(x, p) - E * φ
  x = Roots.find_zero(_to_solve, (0, 1.1*x0), Roots.Bisection(), xrtol=1e-5)                     # first 'rough' estimate
  return Roots.find_zero((_to_solve, xp->der_1_ψ_halo(xp, p)), x, Roots.Newton(),  xrtol=1e-10) # precise estimate with Newton's method
end

# Need to be checked
function phase_space_Eddington(E::Real, x0::Real, p::HaloProfile = nfwProfile; remove_divergence = false) 
  
  # may need to improve the integral computation with methods like double logarithms
  _to_integrate(φ::Real) = der_2_ρ_halo_vs_ψ(x_vs_φ(φ, E, x0, p), p) / sqrt(1-φ)
  integral = QuadGK.quadgk(_to_integrate, 0, 1, rtol=1e-3)
  
  div = 0 # diverging term (non physical)
  !(remove_divergence) && (div = - x0^2 * der_1_ρ_halo(x0, p) / μ_halo(x0, p))
  return 1/(sqrt(8) * π^2) * ( div / sqrt(E) + sqrt(E) * integral )

end 
#-------------


function velocity_dispersion(x::Real, xt::Real, p::αβγProfile = nfwProfile, approx::Bool = true) 
    
    (x >= xt) && (return 0.0)
  
    res = QuadGK.quadgk(lnxp -> ρ_halo(exp(lnxp), p) * μ_halo(exp(lnxp), p) / exp(lnxp), log(x), log(xt), rtol=1e-5)[1]  / ρ_halo(x, p)

    (res < 0) && throw(DomainError("The square of the velocity dispersion cannot be negative (for x = " * string(x) * ", xt = " * string(xt) * ")"))
    
    return sqrt(res)
end

# s- and p-wave luminosity of the halo
λ0_halo(x::Real, p::αβγProfile = nfwProfile) = HypergeometricFunctions._₂F₁((3 - 2*p.γ)/p.α, 2*(p.β - p.γ)/p.α, (3 + p.α - 2*p.γ)/p.α, -x^p.α) * x^(3-2*p.γ) / (3-2*p.γ)
λ1_halo(x::Real, p::αβγProfile = nfwProfile) = -1.0 # TO DO

################################################
#
# Dimensionfull Halo

# Whenever possible do all the computation using the dimensionless profile
# Otherwise unit conversion slightly slow down the code

## Definition of relationships between concentration, mass, scale density and scale radius
function cΔ_from_ρs(ρs::Real, hp::HaloProfile{<:Real}, Δ::Real, ρ_ref::Real)
    g(c::Real) = c^3 / μ_halo(c, hp) - 3 * ρs / Δ / ρ_ref
    return Roots.find_zero(g, (1e-10, 1e+10), Roots.Bisection()) 
end

cΔ_from_ρs(ρs::Real, hp::HaloProfile{<:Real} = nfwProfile, Δ::Real=200, cosmo::Cosmology = plnack18) = cΔ_from_ρs(ρs, hp, Δ, cosmo.bkg.ρ_c0)

mΔ_from_ρs_and_rs(ρs::Real, rs::Real, hp::HaloProfile{<:Real} = nfwProfile, Δ::Real = 200, ρ_ref::Real = planck18_bkg.ρ_c0) = 4 * pi * ρs * rs^3 * μ_halo(cΔ_from_ρs(ρs, hp, Δ, ρ_ref), hp)
ρs_from_cΔ(cΔ::Real, hp::HaloProfile{<:Real} = nfwProfile, Δ::Real = 200, ρ_ref::Real = planck18_bkg.ρ_c0) = Δ * ρ_ref  / 3 * cΔ^3 / μ_halo(cΔ, hp) 
rs_from_cΔ_and_mΔ(cΔ::Real, mΔ::Real, Δ::Real = 200, ρ_ref::Real = planck18_bkg.ρ_c0) =  (3 * mΔ / (4 * π * Δ * ρ_ref))^(1 // 3) / cΔ 
rΔ_from_ρs_and_rs(ρs::Real, rs::Real, hp::HaloProfile{<:Real} = nfwProfile, Δ::Real = 200, ρ_ref::Real = planck18_bkg.ρ_c0) = (3 * mΔ_from_ρs_and_rs(ρs, rs, hp, Δ, ρ_ref) / (4*π*Δ*ρ_ref))^(1//3)

struct Halo{T<:AbstractFloat, S<:Real}
    hp::HaloProfile{S}
    ρs::T
    rs::T
end

Base.length(::Halo) = 1
Base.iterate(iter::Halo) = (iter, nothing)
Base.iterate(::Halo, state::Nothing) = nothing

Base.show(io::IO, h::Halo{<:AbstractFloat, <:Real}) = print(io, "Halo: \n  - " * string(h.hp) * "\n  - ρs = " * string(h.ρs) * " Msun/Mpc^3, rs = " * string(h.rs) * " Mpc \n  - m200 (planck18) = " * string(mΔ(h, 200, planck18)) * " Msun, c200 (planck18) = " * string(cΔ(h, 200, planck18)))

Halo(hp::HaloProfile{S}, ρs::Real, rs::Real) where {S<:Real} = Halo(hp, promote(ρs, rs)...)

function halo_from_mΔ_and_cΔ(hp::HaloProfile{S}, mΔ::T, cΔ::T;  Δ::T = T(200), ρ_ref::T = planck18_bkg.ρ_c0) where {T<:AbstractFloat, S<:Real}
    return Halo{T, S}(hp, convert(T, ρs_from_cΔ(cΔ, hp, Δ, ρ_ref)), convert(T, rs_from_cΔ_and_mΔ(cΔ, mΔ, Δ, ρ_ref)))
end

ρ_halo(r::Real, h::Halo{<:AbstractFloat, <:Real}) = h.ρs * ρ_halo(r/h.rs, h.hp)
m_halo(r::Real, h::Halo{<:AbstractFloat, <:Real}) = 4.0 * π * h.ρs * h.rs^3 * μ_halo(r/h.rs, h.hp)
mΔ(h::Halo{<:AbstractFloat, <:Real}, Δ::Real, ρ_ref::Real) = mΔ_from_ρs_and_rs(h.ρs, h.rs, h.hp, Δ, ρ_ref)
mΔ(h::Halo{<:AbstractFloat, <:Real}, Δ::Real = 200, cosmo::Cosmology = planck18) = mΔ(h, Δ, cosmo.bkg.ρ_c0)
cΔ(h::Halo{<:AbstractFloat, <:Real}, Δ::Real, ρ_ref::Real) = cΔ_from_ρs(h.ρs, h.hp, Δ, ρ_ref)
cΔ(h::Halo{<:AbstractFloat, <:Real}, Δ::Real = 200, cosmo::Cosmology = planck18) = cΔ(h, Δ, cosmo.bkg.ρ_c0)
rΔ(h::Halo{<:AbstractFloat, <:Real}, Δ::Real, ρ_ref::Real) = rΔ_from_ρs_and_rs(h.ρs, h.rs, h.hp, Δ, ρ_ref)
rΔ(h::Halo{<:AbstractFloat, <:Real}, Δ::Real = 200, cosmo::Cosmology = planck18) = rΔ(h, Δ, cosmo.bkg.ρ_c0)

ρs(h::Halo) = h.ρs * ρ_0
rs(h::Halo) = h.rs * r_0

################################################
# Velocity of the particles inside the halo

""" 1D Jeans velocity dispersion in (Mpc / s) """
velocity_dispersion(r::Real, rt::Real, h::Halo) = sqrt(4 * π * h.ρs * h.rs^2 * G_NEWTON) * velocity_dispersion(r/h.rs, rt/h.rs, h.hp) 

""" 1D Jeans velocity dispersion in (km / s) """
velocity_dispersion_kms(r::Real, rt::Real, h::Halo) = velocity_dispersion(r, rt, h) * MPC_TO_KM

""" gravitational potential in (Mpc / s)^2 """
gravitational_potential(r::Real, rt::Real, h::Halo) = - 4 * π * h.ρs * h.rs^2 * G_NEWTON * ψ_halo(r/h.rs, rt/h.rs, h.hp) 

""" gravitational potential in (km / s)^2 """
gravitational_potential_kms(r::Real, rt::Real, h::Halo) = ψ_halo(r, rt, h) * (MPC_TO_KM^2)

""" escape velocity in (Mpc / s) """
escape_velocity(r::Real, rt::Real, h::Halo) = sqrt(2*abs(ψ_halo(r, rt, h)))

""" escape velocity in (km / s) """
escape_velocity_kms(r::Real, rt::Real, h::Halo) = sqrt(2*abs(ψ_halo_kms(r, rt, h))) * MPC_TO_KM

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

halo_from_mΔ(hp::HaloProfile, mΔ::Real, ::Type{T} = SCP12;  Δ::Real = 200, cosmo::Cosmology=planck18, z=0) where {T <: MassConcentrationModel} = halo_from_mΔ_and_cΔ(hp, mΔ, median_concentration(mΔ, T, z = z, cosmo = cosmo), Δ = Δ, ρ_ref = cosmo.bkg.ρ_c0) 

#######################################




function _vd2_nfw(x::Real, xt::Real) 
    (x >= xt) && return 0.0 
    p =(5*x)/3 + (37*x^2)/24 + (11*x^3)/60 - (13*x^4)/240 + (37*x^5)/
  1400 - (17*x^6)/1050 + (991*x^7)/88200 - (1187*x^8)/
  141120 + (15409*x^9)/2328480 - (17929*x^10)/3326400 + (267727*x^11)/
  59459400 - (304027*x^12)/79279200 + 
   (341779*x^13)/103062960 - (190409*x^14)/65585520 + (14314177*x^15)/
  5574769200 - (15715577*x^16)/6861254400 + (108607721*x^17)/
  52766313600 - (39288323*x^18)/21106525440 + (8485835*x^19)/
  5012799792 - (3041361*x^20)/1965803840 + 
   (-(5/3) + (17*x^2)/12 + (11*x^3)/30 - (13*x^4)/120 + (37*x^5)/
     700 - (17*x^6)/525 + (991*x^7)/44100 - (1187*x^8)/
     70560 + (15409*x^9)/1164240 - (17929*x^10)/
     1663200 + (267727*x^11)/29729700 - (304027*x^12)/39639600 + 
        (341779*x^13)/51531480 - (190409*x^14)/
     32792760 + (14314177*x^15)/2787384600 - (15715577*x^16)/
     3430627200 + (108607721*x^17)/26383156800 - (39288323*x^18)/
     10553262720 + (8485835*x^19)/2506399896 - (3041361*x^20)/
     982901920)*xt + 
   (-(37/24) - (17*x)/12 + (11*x^3)/60 - (13*x^4)/240 + (37*x^5)/
     1400 - (17*x^6)/1050 + (991*x^7)/88200 - (1187*x^8)/
     141120 + (15409*x^9)/2328480 - (17929*x^10)/
     3326400 + (267727*x^11)/59459400 - (304027*x^12)/79279200 + 
        (341779*x^13)/103062960 - (190409*x^14)/
     65585520 + (14314177*x^15)/5574769200 - (15715577*x^16)/
     6861254400 + (108607721*x^17)/52766313600 - (39288323*x^18)/
     21106525440 + (8485835*x^19)/5012799792 - (3041361*x^20)/
     1965803840)*xt^2 + 
   (-(11/60) - (11*x)/30 - (11*x^2)/60)*
  xt^3 + (13/240 + (13*x)/120 + (13*x^2)/240)*
  xt^4 + (-(37/1400) - (37*x)/700 - (37*x^2)/1400)*
  xt^5 + (17/1050 + (17*x)/525 + (17*x^2)/1050)*xt^6 + 
   (-(991/88200) - (991*x)/44100 - (991*x^2)/88200)*
  xt^7 + (1187/141120 + (1187*x)/70560 + (1187*x^2)/141120)*
  xt^8 + (-(15409/2328480) - (15409*x)/1164240 - (15409*x^2)/2328480)*
  xt^9 + 
   (17929/3326400 + (17929*x)/1663200 + (17929*x^2)/3326400)*
  xt^10 + (-(267727/59459400) - (267727*x)/29729700 - (267727*x^2)/
     59459400)*
  xt^11 + (304027/79279200 + (304027*x)/39639600 + (304027*x^2)/
     79279200)*xt^12 + 
   (-(341779/103062960) - (341779*x)/51531480 - (341779*x^2)/
     103062960)*
  xt^13 + (190409/65585520 + (190409*x)/32792760 + (190409*x^2)/
     65585520)*
  xt^14 + (-(14314177/5574769200) - (14314177*x)/
     2787384600 - (14314177*x^2)/5574769200)*xt^15 + 
   (15715577/6861254400 + (15715577*x)/3430627200 + (15715577*x^2)/
     6861254400)*
  xt^16 + (-(108607721/52766313600) - (108607721*x)/
     26383156800 - (108607721*x^2)/52766313600)*xt^17 + 
   (39288323/21106525440 + (39288323*x)/10553262720 + (39288323*x^2)/
     21106525440)*
  xt^18 + (-(8485835/5012799792) - (8485835*x)/
     2506399896 - (8485835*x^2)/5012799792)*xt^19 + 
   (3041361/1965803840 + (3041361*x)/982901920 + (3041361*x^2)/
     1965803840)*xt^20
    q = 1/2 + x + xt + 2 * x * xt + 1/2*x^2 + x^2 * xt + 1/2 * xt^2 + x * xt^2 + 1/2*x^2*xt^2
    return x/(1+xt)^2 * ( p + log(xt/x) * q)
end 