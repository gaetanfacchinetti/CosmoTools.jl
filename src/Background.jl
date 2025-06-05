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
#
# Contains functions related to the background cosmology
#
# author: Gaetan Facchinetti
# email: gaetanfacc@gmail.com
#
##################################################################################


export BkgCosmology, FLRW, FlatFLRW, planck18_bkg, edsPlanck18_bkg, convert_cosmo
export hubble_H0, hubble_H0_s, hubble_E, hubble_E2, hubble_H
export Ω, Ω_vs_a, Ω_m, Ω_r, Ω_Λ, mean_ρ, ρ_critical
export z_eq_mr, z_eq_Λm, z_to_a, a_to_z, temperature_CMB_K, k_eq_mr
export growth_factor, growth_factor_Carroll
export lookback_time, lookback_redshift, universe_age
export dflt_bkg_cosmo

@doc raw"""
    BkgCosmology{T<:Real}

Abstract type: generic background cosmology 
"""
abstract type BkgCosmology{T<:AbstractFloat} end

@doc raw"""
    FLRW{T<:Real} <: BkgCosmology{T}

Defines a flat FLRW background cosmology
"""
struct FLRW{T<:AbstractFloat} <: BkgCosmology{T}
    
    # Hubble parameter
    h::T 

    # Abundances of the different components
    Ω_χ0::T
    Ω_b0::T
    Ω_m0::T
    Ω_r0::T
    Ω_γ0::T
    Ω_ν0::T
    Ω_Λ0::T
    Ω_k0::T

    # Derived quantities
    z_eq_mr::T
    z_eq_Λm::T
    k_eq_mr::T
    
    # Quantity with units
    ρ_c0::T
    T0_CMB_K::T

end

Base.length(::BkgCosmology) = 1
Base.iterate(iter::BkgCosmology) = (iter, nothing)
Base.iterate(::BkgCosmology, state::Nothing) = nothing

Base.eltype(::BkgCosmology{T}) where {T<:AbstractFloat} = T

""" Convert bkg_cosmo object attributes to another type """
convert_bkg_cosmo(::Type{T}, bkg_cosmo::FLRW{<:AbstractFloat} = planck18_bkg) where {T<:AbstractFloat} = FLRW{T}([convert(T, getfield(bkg_cosmo, field)) for field in fieldnames(typeof(bkg_cosmo))]...)

# First definition of the abundances
hubble_E2(z::T, Ω_m0::T, Ω_r0::T, Ω_Λ0::T, Ω_k0::T) where {T<:AbstractFloat} = Ω_m0 * (1+z)^3 + Ω_r0 * (1+z)^4 + Ω_Λ0 + Ω_k0 * (1+z)^2

Ω_m(z::T, Ω_m0::T, Ω_r0::T, Ω_Λ0::T, Ω_k0::T) where {T<:AbstractFloat} = Ω_m0 * (1+z)^3 / hubble_E2(z, Ω_m0, Ω_r0, Ω_Λ0, Ω_k0)
Ω_r(z::T, Ω_m0::T, Ω_r0::T, Ω_Λ0::T, Ω_k0::T) where {T<:AbstractFloat} = Ω_r0 * (1+z)^4 / hubble_E2(z, Ω_m0, Ω_r0, Ω_Λ0, Ω_k0)
Ω_Λ(z::T, Ω_m0::T, Ω_r0::T, Ω_Λ0::T, Ω_k0::T) where {T<:AbstractFloat} = Ω_Λ0  / hubble_E2(z, Ω_m0, Ω_r0, Ω_Λ0, Ω_k0)

@doc raw"""
    FLRW(h, Ω_χ0, Ω_b0; T0_CMB_K = 2.72548, Neff = 3.04)

Creates a FLRW instance

# Arguments
- `h::Real` : hubble parameter (dimensionless)
- `Ω_χ0::Real`: cold dark matter abundance today (dimensionless)
- `Ω_b0::Real`: baryon abundance today (dimensionless)
- `T0_CMB_K::Real`: temperature of the CMB today (in Kelvin)
- `Neff::Real`: effective number of neutrinos (dimensionless)

# Returns
- A FLRW structure
"""
function FLRW(h::T, Ω_χ0::T, Ω_b0::T, Ω_k0::T = T(0); T0_CMB_K::T = T(2.72548), Neff::T = T(3.04), EdS::Bool = false) where {T<:AbstractFloat}
    
    # Derived abundances
    Ω_γ0 = EdS ? T(0) : T(4.48131e-7) * T0_CMB_K^4 / h^2
    Ω_ν0 = EdS ? T(0) : Neff * Ω_γ0 * T((7 / 8) * (4 / 11)^(4 / 3))
    Ω_r0 = Ω_γ0 + Ω_ν0
    Ω_m0 = Ω_χ0 + Ω_b0
    Ω_Λ0 = T(1) - Ω_m0 - Ω_r0
 
    ρ_c0 =  3/(8*π*G_NEWTON) * (T(100) * h * KM_TO_MPC)^2

    z_eq_mr = T(0)
    z_eq_Λm = T(0)

    try
        z_eq_mr = EdS ? T(NaN) : exp(Roots.find_zero( y -> Ω_r(exp(y), Ω_m0, Ω_r0, Ω_Λ0, Ω_k0) - Ω_m(exp(y), Ω_m0, Ω_r0, Ω_Λ0, Ω_k0), (T(-10), T(10)), Roots.Bisection())) 
        z_eq_Λm = exp(Roots.find_zero( y -> Ω_Λ(exp(y), Ω_m0, Ω_r0, Ω_Λ0, Ω_k0) - Ω_m(exp(y), Ω_m0, Ω_r0, Ω_Λ0, Ω_k0), (T(-10), T(10)), Roots.Bisection())) 
    catch e
        rethrow("Impossible to definez z_eq_mr and/or z_eq_Λm for this cosmology")
    end

    k_eq_mr =  EdS ? T(NaN) : Ω_r0 / Ω_m0 * (T(100) * h * sqrt(Ω_m0 * (1+z_eq_mr)^3 + Ω_r0 * (1+z_eq_mr)^4 + Ω_Λ0 + Ω_k0 * (1+z_eq_mr)^2) * KM_TO_MPC / C_LIGHT) 

    return FLRW(h, Ω_χ0, Ω_b0, Ω_m0, Ω_r0, Ω_γ0, Ω_ν0, Ω_Λ0, Ω_k0, z_eq_mr, z_eq_Λm, k_eq_mr, ρ_c0, T0_CMB_K) 
end


#########################################################
# Predefinition of some cosmologies
"""
    planck18_bkg::FLRW = FLRW(0.6736, 0.26447, 0.04930)

Flat FLRW cosmology with [*Planck18*](https://arxiv.org/abs/1807.06209) parameters
"""
const planck18_bkg::FLRW{Float64} = FLRW(0.6736, 0.26447, 0.04930)
const planck18_bkg_f32::FLRW{Float32} = convert_bkg_cosmo(Float32, planck18_bkg)
const planck18_bkg_f16::FLRW{Float16} = convert_bkg_cosmo(Float16, planck18_bkg)

"""
    edsplanck18_bkg::FLRW  = FLRW(0.6736, 0.3, 0) 

Flat FLRW cosmology with [*Planck18*](https://arxiv.org/abs/1807.06209) parameters and no baryons
"""
const edsPlanck18_bkg::FLRW{Float64}  = FLRW(0.6736, 0.3, 0.0, EdS = true) 
const edsPlanck18_bkg_f32::FLRW{Float32} = convert_bkg_cosmo(Float32, edsPlanck18_bkg)
const edsPlanck18_bkg_f16::FLRW{Float16} = convert_bkg_cosmo(Float16, edsPlanck18_bkg)


# definition of default cosmology

dflt_bkg_cosmo(::Type{Float64}) = planck18_bkg
dflt_bkg_cosmo(::Type{Float32}) = planck18_bkg_f32
dflt_bkg_cosmo(::Type{Float16}) = planck18_bkg_f16

dflt_bkg_cosmo(::Type{T} = Float64) where {T<:AbstractFloat}  = dflt_bkg_cosmo(T)
#########################################################


#########################################################
# Definition of the hubble parameter and some other related quantities

@doc raw"""
    temperature_CMB_K(z, bkg_cosmo)

CMB temperature (in K) of the Universe at redshift `z`.

# Arguments
- `z::Real`: Redshift (default: 0.0)
- `bkg_cosmo::BkgCosmology`: Cosmological background model (default: `planck18_bkg`)
"""
temperature_CMB_K(z::T, bkg_cosmo::FLRW{T} = dflt_bkg_cosmo(T)) where {T<:AbstractFloat} = bkg_cosmo.T0_CMB_K * (1+z)
temperature_CMB_K(::Type{T} = Float64) where {T<:AbstractFloat} = temperature_CMB_K(T(0))

@doc raw"""
    hubble_H0(bkg_cosmo)

Hubble constant H₀ (in km/s/Mpc) for the cosmology.

# Arguments
- `bkg_cosmo::BkgCosmology`: Cosmological background model (default: `planck18_bkg`)
"""
hubble_H0(bkg_cosmo::FLRW{T}) where {T<:AbstractFloat} = T(100) * bkg_cosmo.h
hubble_H0(::Type{T} = Float64) where {T<:AbstractFloat} = hubble_H0(dflt_bkg_cosmo(T))


@doc raw"""
    hubble_H0_s(bkg_cosmo)

Hubble constant H₀ (in 1/s) for the cosmology.

# Arguments
- `bkg_cosmo::BkgCosmology`: Cosmological background model (default: `planck18_bkg`)
"""
hubble_H0_s(bkg_cosmo::FLRW{T}) where {T<:AbstractFloat} = convert_lengths(hubble_H0(bkg_cosmo), KiloMeters, MegaParsecs)
hubble_H0_s(::Type{T} = Float64) where {T<:AbstractFloat} = hubble_H0_s(dflt_bkg_cosmo(T))

@doc raw"""
    hubble_E2(z, bkg_cosmo)

Squared Hubble evolution function E²(z) (dimensionless) at redshift `z`.

# Arguments
- `z::Real`: Redshift (default: 0.0)
- `bkg_cosmo::BkgCosmology`: Cosmological background model (default: `planck18_bkg`)
"""
hubble_E2(z::T, bkg_cosmo::FLRW{T} = dflt_bkg_cosmo(T)) where {T<:AbstractFloat} = hubble_E2(z, bkg_cosmo.Ω_m0, bkg_cosmo.Ω_r0, bkg_cosmo.Ω_Λ0, bkg_cosmo.Ω_k0)
hubble_E2(::Type{T} = Float64) where {T<:AbstractFloat} = hubble_E2(T(0))


@doc raw"""
    hubble_E(z, bkg_cosmo)

Hubble evolution function E(z) (dimensionless) at redshift `z`.

# Arguments
- `z::Real`: Redshift (default: 0.0)
- `bkg_cosmo::BkgCosmology`: Cosmological background model (default: `planck18_bkg`)
"""
hubble_E(z::T, bkg_cosmo::FLRW{T} = dflt_bkg_cosmo(T)) where {T<:AbstractFloat} = sqrt(hubble_E2(z, bkg_cosmo))
hubble_E(::Type{T} = Float64) where {T<:AbstractFloat} = hubble_E(T(0))


@doc raw"""
    hubble_H(z, bkg_cosmo)

Hubble rate H(z) (in km/s/Mpc) of the Universe at redshift `z`.

# Arguments
- `z::Real`: Redshift (default: 0.0)
- `bkg_cosmo::BkgCosmology`: Cosmological background model (default: `planck18_bkg`)
"""
hubble_H(z::T, bkg_cosmo::FLRW{T} = dflt_bkg_cosmo(T)) where {T<:AbstractFloat} = hubble_E(z, bkg_cosmo) .* hubble_H0(bkg_cosmo)
hubble_H(::Type{T} = Float64) where {T<:AbstractFloat} = hubble_H(T(0))


@doc raw"""
    k_eq_mr(bkg_cosmo)

Wave number at matter–radiation equality.

# Arguments
- `bkg_cosmo::BkgCosmology`: Cosmological background model
"""
k_eq_mr(bkg_cosmo::FLRW) = bkg_cosmo.k_eq_mr
#########################################################


#########################################################
# Definition of the densities and abundances
# All densities are in units of Msun / Mpc^3 

export Species, Radiation, Photons, Neutrinos, Matter, ColdDarkMatter, Baryons, DarkEnergy, Curvature

@doc raw"""
    Species

Abstract type: generic species populating the Universe
"""
abstract type Species end

@doc raw"""
    Radiation <: Species
"""
abstract type Radiation <: Species end

@doc raw"""
    Photons <: Radiation
"""
abstract type Photons <: Radiation end

@doc raw"""
    Neutrinos <: Radiation
"""
abstract type Neutrinos <: Radiation end

@doc raw"""
    Matter <: Species
"""
abstract type Matter <: Species end

@doc raw"""
    ColdDarkMatter <: Matter
"""
abstract type ColdDarkMatter <: Matter end

@doc raw"""
    Baryons <: Matter
"""
abstract type Baryons <: Matter end

@doc raw"""
    DarkEnergy <: Species
"""
abstract type DarkEnergy <: Species end

@doc raw"""
    Curvature <: Species
"""
abstract type Curvature <: Species end


Ω(z::T, n::Int, Ω0::T, bkg_cosmo::FLRW{T}) where {T<:AbstractFloat} = Ω0 * (1+z)^n / hubble_E2(z, bkg_cosmo)

Ω0(::Type{Matter}, bkg_cosmo::FLRW)         = bkg_cosmo.Ω_m0
Ω0(::Type{Radiation}, bkg_cosmo::FLRW)      = bkg_cosmo.Ω_r0
Ω0(::Type{DarkEnergy}, bkg_cosmo::FLRW)     = bkg_cosmo.Ω_Λ0
Ω0(::Type{Photons}, bkg_cosmo::FLRW)        = bkg_cosmo.Ω_γ0
Ω0(::Type{Neutrinos}, bkg_cosmo::FLRW)      = bkg_cosmo.Ω_ν0
Ω0(::Type{ColdDarkMatter}, bkg_cosmo::FLRW) = bkg_cosmo.Ω_χ0
Ω0(::Type{Baryons}, bkg_cosmo::FLRW)        = bkg_cosmo.Ω_b0
Ω0(::Type{Curvature}, bkg_cosmo::FLRW)      = bkg_cosmo.Ω_k0

index_ρ(::Type{S}) where {S<:Matter} = 3
index_ρ(::Type{S}) where {S<:Radiation}   = 4
index_ρ(::Type{DarkEnergy})::Int  = 0

@doc raw"""
    Ω(S, z = 0, bkg_cosmo = planck18_bkg) where {S<:Species}

Abundance of species `S` at redshift `z` and for the background cosmology `bkg_cosmo`
"""
Ω(::Type{S}, z::T, bkg_cosmo::FLRW{T} = dflt_bkg_cosmo(T)) where {S<:Species, T<:AbstractFloat} = Ω(z, index_ρ(S), Ω0(S, bkg_cosmo), bkg_cosmo)
Ω(::Type{S}, ::Type{T} = Float64) where {S<:Species, T<:AbstractFloat} = Ω(S, T(0))

Ω_vs_a(::Type{S}, a::T, bkg_cosmo::FLRW{T} = dflt_bkg_cosmo(T)) where {S<:Species, T<:AbstractFloat} = Ω(S, a_to_z(a), bkg_cosmo)


""" Critical density (in Msun/Mpc^3) of the Universe at redshift `z` (by default z=0) for the cosmology `bkg_cosmo` """
ρ_critical(z::T, bkg_cosmo::FLRW{T} = dflt_bkg_cosmo(T)) where {T<:AbstractFloat} = bkg_cosmo.ρ_c0 * hubble_E2(z, bkg_cosmo)
ρ_critical(::Type{T} = Float64) where {T<:AbstractFloat} = ρ_critical(T(0))

mean_ρ(::Type{Radiation}, z::T, bkg_cosmo::FLRW{T}) where {T<:AbstractFloat}      = bkg_cosmo.Ω_r0 * bkg_cosmo.ρ_c0 * (1+z)^4
mean_ρ(::Type{Photons}, z::T, bkg_cosmo::FLRW{T}) where {T<:AbstractFloat}        = bkg_cosmo.Ω_γ0 * bkg_cosmo.ρ_c0 * (1+z)^4
mean_ρ(::Type{Neutrinos}, z::T, bkg_cosmo::FLRW{T}) where {T<:AbstractFloat}      = bkg_cosmo.Ω_ν0 * bkg_cosmo.ρ_c0 * (1+z)^4
mean_ρ(::Type{ColdDarkMatter}, z::T, bkg_cosmo::FLRW{T}) where {T<:AbstractFloat} = bkg_cosmo.Ω_χ0 * bkg_cosmo.ρ_c0 * (1+z)^3
mean_ρ(::Type{Matter}, z::T, bkg_cosmo::FLRW{T}) where {T<:AbstractFloat}         = bkg_cosmo.Ω_m0 * bkg_cosmo.ρ_c0 * (1+z)^3
mean_ρ(::Type{Baryons}, z::T, bkg_cosmo::FLRW{T}) where {T<:AbstractFloat}        = bkg_cosmo.Ω_b0 * bkg_cosmo.ρ_c0 * (1+z)^3
mean_ρ(::Type{Curvature}, z::T, bkg_cosmo::FLRW{T}) where {T<:AbstractFloat}      = bkg_cosmo.Ω_k0 * bkg_cosmo.ρ_c0 * (1+z)^2
mean_ρ(::Type{DarkEnergy}, z::T, bkg_cosmo::FLRW{T}) where {T<:AbstractFloat}     = bkg_cosmo.Ω_Λ0 * bkg_cosmo.ρ_c0

@doc raw"""
    mean_ρ(T, z = 0, bkg_cosmo = planck18_bkg) where {T<:Species}

Energy density of species `T` at redshift `z` and for the background cosmology `bkg_cosmo`
"""
mean_ρ(::Type{S}, z::T, bkg_cosmo::FLRW{T} = dflt_bkg_cosmo(T)) where {S<:Species, T<:AbstractFloat} = mean_ρ(S, z, bkg_cosmo)
mean_ρ(::Type{S}, ::Type{T}) where {S<:Species, T<:AbstractFloat} = mean_ρ(S, T(0), bkg_cosmo)


#########################################################


#########################################################
# Definition of specific time quantities

@doc raw""" 
    z_to_a(z)

Convert redshift z to scale parameter a. Simple implementation of

``a = \frac{1}{1+z}``

# Arguments
- `z::Real` : redshift (dimensionless)
"""
z_to_a(z::AbstractFloat) = 1 / (1 + z)


@doc raw""" 
    a_to_z(z)

Convert scale parameter to redshift. Simple implementation of

``z = \frac{1}{a}-1``

# Arguments
- `a::Real` : scale factor (dimensionless)
"""
a_to_z(a::AbstractFloat) = 1 / a - 1

z_eq_mr(bkg_cosmo::FLRW) = bkg_cosmo.z_eq_mr
z_eq_Λm(bkg_cosmo::FLRW) = bkg_cosmo.z_eq_Λm
a_eq_mr(bkg_cosmo::FLRW) = z_to_a(bkg_cosmo.z_eq_mr)
a_eq_Λm(bkg_cosmo::FLRW) = z_to_a(bkg_cosmo.z_eq_Λm)

z_eq_mr(::Type{T} = Float64) where {T<:AbstractFloat} = z_eq_mr(dflt_bkg_cosmo(T))
z_eq_Λm(::Type{T} = Float64) where {T<:AbstractFloat} = z_eq_Λm(dflt_bkg_cosmo(T))
a_eq_mr(::Type{T} = Float64) where {T<:AbstractFloat} = a_eq_mr(dflt_bkg_cosmo(T))
a_eq_Λm(::Type{T} = Float64) where {T<:AbstractFloat} = a_eq_Λm(dflt_bkg_cosmo(T))

δt_s(a0::T, a1::T, bkg_cosmo::FLRW{T} = dflt_bkg_cosmo(T); kws...) where {T<:AbstractFloat} = QuadGK.quadgk(a -> 1 / hubble_E(a_to_z(a), bkg_cosmo) / a, a0, a1, rtol=T(1e-3); kws...)[1] / (hubble_H0_s(bkg_cosmo))

@doc raw""" 

    universe_age(z, bkg_cosmo; kws...)
    
age of the universe (s)

# Arguments
- `z::Float64` : redshift (dimensionless), default 0.0
- `bkg_cosmo::BkgCosmology` : background cosmology (dimensionless), default planck18_bkg

# Kwargs
- Any argument that can be passed to `QuadGK.quadgk'
"""
universe_age(z::T, bkg_cosmo::FLRW{T} = dflt_bkg_cosmo(T); kws...) where {T<:AbstractFloat} = δt_s(T(0), z_to_a(z), bkg_cosmo; kws...)
universe_age(::Type{T} = Float64; kws...) where {T<:AbstractFloat} = universe_age(T(0); kws...)

@doc raw""" 

    lookback_time(z, bkg_cosmo; kws...)
    
lookback time of the Universe (s)

# Arguments
- `z::Real` : redshift (dimensionless)
- `bkg_cosmo::BkgCosmology` : background cosmology (dimensionless), default planck18_bkg

# Kwargs
- Any argument that can be passed to `QuadGK.quadgk'
"""
lookback_time(z::T, bkg_cosmo::FLRW{T} = dflt_bkg_cosmo(T); kws...) where {T<:AbstractFloat} = δt_s(z_to_a(z), T(1), bkg_cosmo; kws...)

@doc raw""" 

    lookback_lookback_redshift(t, bkg_cosmo; kws...)
    
lookback redshift of the Universe for t in (s)

# Arguments
- `t::Real` : time (s)
- `bkg_cosmo::BkgCosmology` : background cosmology (dimensionless), default planck18_bkg

# Kwargs
- Any argument that can be passed to `QuadGK.quadgk'
"""
lookback_redshift(t::T, bkg_cosmo::FLRW{T} = dflt_bkg_cosmo(T); kws...) where {T<:AbstractFloat} = exp10(Roots.find_zero(log10_z -> lookback_time(exp10(log10_z), bkg_cosmo; kws...) - t, (T(-15), T(15)), Roots.Bisection())) 
#########################################################


#########################################################
# Functions related to perturbations
"""
    growth_factor(z, bkg_cosmo)

    Exact growth factor in a matter-Λ Universe computed from an integral
    Carroll et al. 1992 (Mo and White p. 172)
    Corresponds to D1(a=1) with the definition of Dodelson 2003

# Arguments
- `z::Real` : redshift (dimensionless)
- `bkg_cosmo::BkgCosmology` : background cosmology, default planck18_bkg

# Kwargs
- Any argument that can be passed to `QuadGK.quadgk'
"""
function growth_factor(z::T, bkg_cosmo::FLRW{T} = dflt_bkg_cosmo(T); rtol=T(1e-6), kws...) where {T<:AbstractFloat}
    norm = T(2.5) * bkg_cosmo.Ω_m0 * sqrt(bkg_cosmo.Ω_m0 * (1+z)^3 + bkg_cosmo.Ω_Λ0)
    return norm * QuadGK.quadgk(a -> (bkg_cosmo.Ω_m0 * a^(-1) + bkg_cosmo.Ω_Λ0 * a^2)^(-3/2), T(0), z_to_a(z), rtol=rtol; kws...)[1]
end

"""
    growth_factor_Carroll(z, bkg_cosmo)

    Approximate growth factor in a matter-Λ Universe (faster than growth_factor)
    Carroll et al. 1992 (Mo and White p. 172)
    Corresponds to D1(a=1) with the definition of Dodelson 2003

# Arguments
- `z::Real` : redshift (dimensionless)
- `bkg_cosmo::BkgCosmology` : background cosmology, default planck18_bkg
"""
function growth_factor_Carroll(z::T, bkg_cosmo::FLRW{T} = dflt_bkg_cosmo(T)) where {T<:AbstractFloat}
    # Abundances in a Universe without radiation
    _Ω_m = bkg_cosmo.Ω_m0 * (1+z)^3 / (bkg_cosmo.Ω_m0 * (1+z)^3 + bkg_cosmo.Ω_Λ0)
    _Ω_Λ = bkg_cosmo.Ω_Λ0 / (bkg_cosmo.Ω_m0 * (1+z)^3 + bkg_cosmo.Ω_Λ0)
    return T(2.5)*_Ω_m/(_Ω_m^(T(4.0/7.0)) - _Ω_Λ + (1 + _Ω_m/2) * (1 + _Ω_Λ/70))/(1+z)
end
#########################################################


