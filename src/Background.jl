########################################################
# Subfile of the module Cosmojuly.jl
#
# Contains functions related to the transfer functions
#
# author: Gaetan Facchinetti
# email: gaetanfacc@gmail.com
#
########################################################


# First define some basic constants (the only one we will need)
const G_NEWTON::Float64 = 4.517103049894965e-48 # in Mpc^3 / Msun / s^(2)
const C_LIGHT::Float64 = 9.715611890180198e-15 # Mpc / s 
const MPC_TO_KM::Float64 = 3.085677581491367e19
const KM_TO_MPC::Float64 = 1.0/MPC_TO_KM

export BkgCosmology, FLRW, FlatFLRW, planck18_bkg, edsPlanck18_bkg, convert_cosmo
export hubble_H0, hubble_E, hubble_E2
export Ω, Ω_vs_a, Ω_m, Ω_r, Ω_Λ, mean_ρ, ρ_critical
export z_eq_mr, z_eq_Λm, z_to_a, a_to_z, temperature_CMB_K, k_eq_mr
export growth_factor, growth_factor_Carroll
export lookback_time, lookback_redshift, universe_age

@doc raw"""
    BkgCosmology

Abstract type: generic background cosmology 
"""
abstract type BkgCosmology end

@doc raw"""
    FLRW{T<:Real} <: BkgCosmology

Defines a flat FLRW background cosmology
"""
struct FLRW{T<:Real} <: BkgCosmology
    
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

""" Convert bkg_cosmo object attributes to another type """
convert_cosmo(::Type{T}, bkg_cosmo::FLRW{<:Real} = planck18_bkg) where {T<:Real} = FLRW{T}([convert(T, getfield(bkg_cosmo, field)) for field in fieldnames(typeof(bkg_cosmo))]...)

# First definition of the abundances
hubble_E2(z::Real, Ω_m0::Real, Ω_r0::Real, Ω_Λ0::Real, Ω_k0::Real) = Ω_m0 * (1+z)^3 + Ω_r0 * (1+z)^4 + Ω_Λ0 + Ω_k0 * (1+z)^2

Ω_m(z::Real, Ω_m0::Real, Ω_r0::Real, Ω_Λ0::Real, Ω_k0::Real) = Ω_m0 * (1+z)^3 / hubble_E2(z, Ω_m0, Ω_r0, Ω_Λ0, Ω_k0)
Ω_r(z::Real, Ω_m0::Real, Ω_r0::Real, Ω_Λ0::Real, Ω_k0::Real) = Ω_r0 * (1+z)^4 / hubble_E2(z, Ω_m0, Ω_r0, Ω_Λ0, Ω_k0)
Ω_Λ(z::Real, Ω_m0::Real, Ω_r0::Real, Ω_Λ0::Real, Ω_k0::Real) = Ω_Λ0  / hubble_E2(z, Ω_m0, Ω_r0, Ω_Λ0, Ω_k0)

@doc raw"""
    FlatFLRW(h, Ω_χ0, Ω_b0; T0_CMB_K = 2.72548, Neff = 3.04)

Creates a FlatFLRW instance

# Arguments
- `h::Real` : hubble parameter (dimensionless)
- `Ω_χ0::Real`: cold dark matter abundance today (dimensionless)
- `Ω_b0::Real`: baryon abundance today (dimensionless)
- `T0_CMB_K::Real`: temperature of the CMB today (in Kelvin)
- `Neff::Real`: effective number of neutrinos (dimensionless)
"""
function FLRW(h::Real, Ω_χ0::Real, Ω_b0::Real, Ω_k0::Real=0; T0_CMB_K::Real = 2.72548, Neff::Real = 3.04, EdS::Bool = false)
    
    # Derived abundances
    Ω_γ0 = EdS ? 0.0 : 4.48131e-7 * T0_CMB_K^4 / h^2
    Ω_ν0 = EdS ? 0.0 : Neff * Ω_γ0 * (7 / 8) * (4 / 11)^(4 / 3)
    Ω_r0 = Ω_γ0 + Ω_ν0
    Ω_m0 = Ω_χ0 + Ω_b0
    Ω_Λ0 = 1.0 - Ω_m0 - Ω_r0
 
    ρ_c0 =  3/(8*π*G_NEWTON) * (100 * h * KM_TO_MPC)^2

    z_eq_mr = 0.0
    z_eq_Λm = 0.0

    try
        z_eq_mr = EdS ? NaN : exp(Roots.find_zero( y -> Ω_r(exp(y), Ω_m0, Ω_r0, Ω_Λ0, Ω_k0) - Ω_m(exp(y), Ω_m0, Ω_r0, Ω_Λ0, Ω_k0), (-10, 10), Roots.Bisection())) 
        z_eq_Λm = exp(Roots.find_zero( y -> Ω_Λ(exp(y), Ω_m0, Ω_r0, Ω_Λ0, Ω_k0) - Ω_m(exp(y), Ω_m0, Ω_r0, Ω_Λ0, Ω_k0), (-10, 10), Roots.Bisection())) 
    catch e
        Throw(e("Impossible to definez z_eq_mr and/or z_eq_Λm for this cosmology"))

    end

    k_eq_mr =  EdS ? NaN : Ω_r0 / Ω_m0 * (100. * h * sqrt(Ω_m0 * (1+z_eq_mr)^3 + Ω_r0 * (1+z_eq_mr)^4 + Ω_Λ0 + Ω_k0 * (1+z_eq_mr)^2) * KM_TO_MPC / C_LIGHT) 

    return FLRW(promote(h, Ω_χ0, Ω_b0, Ω_m0, Ω_r0, Ω_γ0, Ω_ν0, Ω_Λ0, Ω_k0, z_eq_mr, z_eq_Λm, k_eq_mr, ρ_c0, T0_CMB_K)...) 
end


#########################################################
# Predefinition of some cosmologies
"""
    planck18_bkg::FlatFLRW = FlatFLRW(0.6736, 0.26447, 0.04930)

Flat FLRW cosmology with [*Planck18*](https://arxiv.org/abs/1807.06209) parameters
"""
const planck18_bkg::FLRW = FLRW(0.6736, 0.26447, 0.04930)

"""
    edsplanck18_bkg::FlatFLRW  = FlatFLRW(0.6736, 0.3, 0) 

Flat FLRW cosmology with [*Planck18*](https://arxiv.org/abs/1807.06209) parameters and no baryons
"""
const edsPlanck18_bkg::FLRW  = FLRW(0.6736, 0.3, 0, EdS = true) 
#########################################################


#########################################################
# Definition of the hubble parameter and some other related quantities

""" CMB temperature (in K) of the Universe at redshift `z` (by default z=0) for the cosmology `bkg_cosmo` """
temperature_CMB_K(z::Real = 0.0, bkg_cosmo::BkgCosmology = planck18_bkg)  = bkg_cosmo.T0_CMB_K * (1+z)

""" Hubble constant H0 (in km/s/Mpc) for the cosmology `bkg_cosmo` """
hubble_H0(bkg_cosmo::BkgCosmology = planck18_bkg) = 100 * bkg_cosmo.h

""" Hubble constant H0 (in 1/s) for the cosmology `bkg_cosmo` """
hubble_H0_s(bkg_cosmo::BkgCosmology = planck18_bkg) = hubble_H0(bkg_cosmo) * KM_TO_MPC

""" Hubble evolution (no dimension) of the Universe squared at redshift `z` (by default z=0) for the cosmology `bkg_cosmo` """
hubble_E2(z::Real = 0.0, bkg_cosmo::BkgCosmology = planck18_bkg) = hubble_E2(z, bkg_cosmo.Ω_m0, bkg_cosmo.Ω_r0, bkg_cosmo.Ω_Λ0, bkg_cosmo.Ω_k0)

""" Hubble evolution (no dimension) of the Universe at redshift `z` (by default z=0) for the cosmology `bkg_cosmo` """
hubble_E(z::Real = 0.0, bkg_cosmo::BkgCosmology = planck18_bkg) = sqrt(hubble_E2(z, bkg_cosmo))

""" Hubble rate H(z) (in km/s/Mpc) of the Universe at redshift `z` (by default z=0) for the cosmology `bkg_cosmo` """
hubble_H(z::Real = 0.0, bkg_cosmo::BkgCosmology = planck18_bkg) = hubble_E(z, bkg_cosmo) .* hubble_H0(bkg_cosmo) 

""" k at equivalence between matter and radiation """
k_eq_mr(bkg_cosmo::BkgCosmology) = bkg_cosmo.k_eq_mr
#########################################################


#########################################################
# Definition of the densities and abundances
# All densities are in units of Msun / Mpc^3 

export Species, Radiation, Photons, Neutrinos, Matter, ColdDarkMatter, Baryons, DarkEnergy, Curvature

abstract type Species end
abstract type Radiation <: Species end
abstract type Photons <: Radiation end
abstract type Neutrinos <: Radiation end
abstract type Matter <: Species end
abstract type ColdDarkMatter <: Matter end
abstract type Baryons <: Matter end
abstract type DarkEnergy <: Species end
abstract type Curvature <: Species end

Ω(z::Real, n::Int, Ω0::Real, bkg_cosmo::BkgCosmology) = Ω0 * (1+z)^n / hubble_E2(z, bkg_cosmo)

Ω0(::Type{Matter}, bkg_cosmo::BkgCosmology)         = bkg_cosmo.Ω_m0
Ω0(::Type{Radiation}, bkg_cosmo::BkgCosmology)      = bkg_cosmo.Ω_r0
Ω0(::Type{DarkEnergy}, bkg_cosmo::BkgCosmology)     = bkg_cosmo.Ω_Λ0
Ω0(::Type{Photons}, bkg_cosmo::BkgCosmology)        = bkg_cosmo.Ω_γ0
Ω0(::Type{Neutrinos}, bkg_cosmo::BkgCosmology)      = bkg_cosmo.Ω_ν0
Ω0(::Type{ColdDarkMatter}, bkg_cosmo::BkgCosmology) = bkg_cosmo.Ω_χ0
Ω0(::Type{Baryons}, bkg_cosmo::BkgCosmology)        = bkg_cosmo.Ω_b0
Ω0(::Type{Curvature}, bkg_cosmo::BkgCosmology)      = bkg_cosmo.Ω_k0

index_ρ(::Type{T}) where {T<:Matter} = 3
index_ρ(::Type{T}) where {T<:Radiation}   = 4
index_ρ(::Type{DarkEnergy})::Int  = 0

@doc raw"""
    Ω(T, z = 0, bkg_cosmo = planck18_bkg) where {T<:Species}

Abundance of species `T` at redshift `z` and for the background cosmology `bkg_cosmo`
"""
Ω(::Type{T}, z::Real = 0, bkg_cosmo::BkgCosmology = planck18_bkg) where {T<:Species} = Ω(z, index_ρ(T), Ω0(T, bkg_cosmo), bkg_cosmo)
Ω_vs_a(::Type{T}, a::Real = 1, bkg_cosmo::BkgCosmology = planck18_bkg) where {T<:Species} = Ω(T, a_to_z(a), bkg_cosmo)


""" Critical density (in Msun/Mpc^3) of the Universe at redshift `z` (by default z=0) for the cosmology `bkg_cosmo` """
ρ_critical(z::Real = 0, bkg_cosmo::BkgCosmology = planck18_bkg) = bkg_cosmo.ρ_c0 * hubble_E2(z, bkg_cosmo)

mean_ρ(::Type{Radiation}, z::Real, bkg_cosmo::BkgCosmology)      = bkg_cosmo.Ω_r0 * bkg_cosmo.ρ_c0 * (1+z)^4
mean_ρ(::Type{Photons}, z::Real, bkg_cosmo::BkgCosmology)        = bkg_cosmo.Ω_γ0 * bkg_cosmo.ρ_c0 * (1+z)^4
mean_ρ(::Type{Neutrinos}, z::Real, bkg_cosmo::BkgCosmology)      = bkg_cosmo.Ω_ν0 * bkg_cosmo.ρ_c0 * (1+z)^4
mean_ρ(::Type{ColdDarkMatter}, z::Real, bkg_cosmo::BkgCosmology) = bkg_cosmo.Ω_χ0 * bkg_cosmo.ρ_c0 * (1+z)^3
mean_ρ(::Type{Matter}, z::Real, bkg_cosmo::BkgCosmology)         = bkg_cosmo.Ω_m0 * bkg_cosmo.ρ_c0 * (1+z)^3
mean_ρ(::Type{Baryons}, z::Real, bkg_cosmo::BkgCosmology)        = bkg_cosmo.Ω_b0 * bkg_cosmo.ρ_c0 * (1+z)^3
mean_ρ(::Type{Curvature}, z::Real, bkg_cosmo::BkgCosmology)      = bkg_cosmo.Ω_k0 * bkg_cosmo.ρ_c0 * (1+z)^2
mean_ρ(::Type{DarkEnergy}, z::Real, bkg_cosmo::BkgCosmology)     = bkg_cosmo.Ω_Λ0 * bkg_cosmo.ρ_c0

@doc raw"""
    mean_ρ(T, z = 0, bkg_cosmo = planck18_bkg) where {T<:Species}

Energy density of species `T` at redshift `z` and for the background cosmology `bkg_cosmo`
"""
mean_ρ(::Type{T}, z::Real = 0, bkg_cosmo::BkgCosmology = planck18_bkg) where {T<:Species} = mean_ρ(T, z, bkg_cosmo)

#########################################################


#########################################################
# Definition of specific time quantities

""" Convert redshift to scale parameter """
z_to_a(z::Real) = 1 / (1 + z)
""" Convert scale parameter to redshift """
a_to_z(a::Real) = 1 / a - 1

z_eq_mr(bkg_cosmo::BkgCosmology = planck18_bkg) = bkg_cosmo.z_eq_mr
z_eq_Λm(bkg_cosmo::BkgCosmology = planck18_bkg) = bkg_cosmo.z_eq_Λm
a_eq_mr(bkg_cosmo::BkgCosmology = planck18_bkg) = z_to_a(bkg_cosmo.z_eq_mr)
a_eq_Λm(bkg_cosmo::BkgCosmology = planck18_bkg) = z_to_a(bkg_cosmo.z_eq_Λm)

δt_s(a0::Real, a1::Real, bkg_cosmo::BkgCosmology = planck18_bkg; kws...) = QuadGK.quadgk(a -> 1 / hubble_E(a_to_z(a), bkg_cosmo) / a, a0, a1, rtol=1e-3; kws...)[1] / (hubble_H0(bkg_cosmo) * km_to_Mpc)

""" age of the universe (s)"""
universe_age(z::Float64 = 0.0, bkg_cosmo::BkgCosmology = planck18_bkg; kws...) = δt_s(0, z_to_a(z), bkg_cosmo; kws...)

""" lookback time of the Universe (s) """
lookback_time(z::Real, bkg_cosmo::BkgCosmology = planck18_bkg; kws...) = δt_s(z_to_a(z), 1, bkg_cosmo; kws...)

""" lookback redshift of the Universe for t in (s) """
lookback_redshift(t::Real, bkg_cosmo::BkgCosmology = planck18_bkg; kws...) = exp(Roots.find_zero(lnz -> lookback_time(exp(lnz), bkg_cosmo; kws...)-t, (-5, 10), Roots.Bisection())) 
#########################################################


#########################################################
# Functions related to perturbations
"""
    growth_factor(z, bkg_cosmo)

    Exact growth factor in a matter-Λ Universe computed from an integral
    Carroll et al. 1992 (Mo and White p. 172)
    Corresponds to D1(a=1) with the definition of Dodelson 2003

# Arguments
- z: redshift
- bkg_cosmo: background cosmology (default Planck18)
"""
function growth_factor(z::Real, bkg_cosmo::BkgCosmology = planck18_bkg; rtol=1e-6, kws...)::Real
    norm = 2.5 * bkg_cosmo.Ω_m0 * sqrt(bkg_cosmo.Ω_m0 * (1+z)^3 + bkg_cosmo.Ω_Λ0)
    return norm * QuadGK.quadgk(a -> (bkg_cosmo.Ω_m0 * a^(-1) + bkg_cosmo.Ω_Λ0 * a^2)^(-3/2), 0, z_to_a(z), rtol=rtol; kws...)[1]
end

"""
    growth_factor_Carroll(z, bkg_cosmo)

    Approximate growth factor in a matter-Λ Universe (faster than growth_factor)
    Carroll et al. 1992 (Mo and White p. 172)
    Corresponds to D1(a=1) with the definition of Dodelson 2003

# Arguments
- z: redshift
- bkg_cosmo: background cosmology (default Planck18)
"""
function growth_factor_Carroll(z::Real, bkg_cosmo=BkgCosmology = planck18_bkg)::Real
    # Abundances in a Universe without radiation
    _Ω_m = bkg_cosmo.Ω_m0 * (1+z)^3 / (bkg_cosmo.Ω_m0 * (1+z)^3 + bkg_cosmo.Ω_Λ0)
    _Ω_Λ = bkg_cosmo.Ω_Λ0 / (bkg_cosmo.Ω_m0 * (1+z)^3 + bkg_cosmo.Ω_Λ0)
    return 2.5*_Ω_m/(_Ω_m^(4.0/7.0) - _Ω_Λ + (1.0 + 0.5*_Ω_m) * (1.0 + 1.0/70.0*_Ω_Λ))/(1+z)
end
#########################################################


