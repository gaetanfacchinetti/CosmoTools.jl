export Cosmology, planck18, curvature_power_spectrum, matter_power_spectrum, power_spectrum_ΛCDM
export window_function, Window, TopHatWindow, SharpKWindow, GaussianWindow, radius_from_mass, mass_from_radius, dradius_dmass
export σ2_mps, dσ2_mps_dR, σ_mps, dσ_mps_dR, σ2_mps_M, dσ2_mps_dM, σ_mps_M, dσ_mps_dM

#####################
# COSMOLOGY STRUCTURE
#####################


@doc raw"""
   Cosmology

Defines a generic cosmology
"""
struct Cosmology
    name::String
    bkg::BkgCosmology
    power_spectrum::Function
    transfer_function_model::TransferFunctionModel
end


function Cosmology(name::String, bkg::BkgCosmology, power_spectrum::Function, ::Type{T} = EH98) where {T<:TransferFunctionModel}
    return Cosmology(name, bkg, power_spectrum, T(bkg))
end

const planck18::Cosmology = Cosmology("PLANCK18", planck18_bkg, k->power_spectrum_ΛCDM(k, 1e-10*exp(3.044), 0.9649), EH98_planck18)


######################################
# (CURVATURE AND MATTER) POWER SPECTRA
######################################

""" ΛCDM power-law power spectrum (dimensionless) at k (in 1/Mpc) """
function power_spectrum_ΛCDM(k::Real, amplitude::Real = 1e-10*exp(3.044), index::Real = 0.9649)
    return amplitude * (k / 0.05)^(index-1)
end

""" Curvature power spectrum (in Mpc^3) at k (in 1/Mpc) """
function curvature_power_spectrum(k::Real, power_spectrum::Function = power_spectrum_ΛCDM)
    return 2.0 * pi^2 * power_spectrum(k) / k^3
end 


""" Matter power spectrum (in Mpc^3) at k (in 1/Mpc) """
function matter_power_spectrum( 
    k::Real, 
    z::Real = 0.0;
    cosmology::Cosmology = planck18,
    growth_function::Function = growth_factor_Carroll,
    dimensionless = false,
    with_baryons::Bool = true)

    _c_over_H0_Mpc = C_LIGHT / hubble_H0_s(cosmology.bkg)
    _D1_z = growth_function(z, cosmology.bkg)
    _tf = transfer_function(k, cosmology.transfer_function_model, with_baryons = with_baryons)
    _prefactor = !dimensionless ? 1 : k^3/(2*π^2)

    return _prefactor * (4. / 25.) * (_D1_z * k^2 * _tf * _c_over_H0_Mpc^2 / cosmology.bkg.Ω_m0)^2 * curvature_power_spectrum(k, cosmology.power_spectrum) 
end


###################################################
# QUANTITIES DERIVED FROM THE MATTER POWER SPECTRUM
###################################################

abstract type Window end
abstract type TopHatWindow <: Window end
abstract type SharpKWindow <: Window end
abstract type GaussianWindow <: Window end

# -------------------------------
# WINDOW FUNCTION AND DERIVATIVES

@doc raw""" 
    window_function(kR, [T])

Give the window function for the product of mode by radius ``k \times R``: `kR::Real` and a certain window type `T<:Window` 
"""
window_function(kR::Real, ::Type{T} = TopHatWindow) where {T<:Window} = window_function(kR, T)

window_function(kR::Real, ::Type{TopHatWindow})      = kR > 1e-3 ? 3.0 * (sin(kR) - kR * cos(kR)) / kR^3 : 1.0 - kR^2/10.0 
window_function(kR::Real, ::Type{SharpKWindow})      = 1 - kR > 0 ? 1.0 : 0.0
window_function(kR::Real, ::Type{GaussianWindow}) = exp(-kR * kR / 2.0)

dwindow_function_dkR(kR::Real, ::Type{TopHatWindow})   = kR < 0.1 ? -kR/5.0 + kR^3/70.0 : 3*sin(kR) / kR^2 - 9 * (-kR * cos(kR) + sin(kR)) / kR^4
dwindow_function_dkR(kR::Real,  ::Type{SharpKWindow})   = 0
dwindow_function_dkR(kR::Real, ::Type{GaussianWindow}) = -kR * exp(-kR^2/2.0)

# --------------
# VOLUME FACTORS

volume_factor(::Type{TopHatWindow})      = 4.0/3.0 * π 
volume_factor(::Type{SharpKWindow})      = 6.0 * π^2 
volume_factor(::Type{GaussianWindow}) = (2*π)^(3//2)

@doc raw""" 
   volume_factor([T])

Give the volume factor associated to a certain window type `T<:Window` (default is `TopHatWindow`)
"""
volume_factor(::Type{T} = TopHatWindow) where {T<:Window} = volume_factor(T)

# ------------------------------------
# MASS - LAGRANGIAN RADIUS CONVERSIONS

""" 
    mass_from_radius(r, [T, [cosmo]]) in Msun

Give the Lagrangian mass (in Msun) in terms of the comoving radius R (in Mpc)

# Arguments
- `r`: radius in Mpc
- `T`: type of Window (default is `TopHatWindow`)
- `cosmo` : cosmology type (default is `planck18_bkg`)
"""
mass_from_radius(r::Real, ::Type{T} = TopHatWindow, cosmo::BkgCosmology = planck18_bkg) where {T<:Window}  = volume_factor(T) * mean_ρ(Matter, 0, cosmo) * r^3


@doc raw""" 
    radius_from_mass(m, [T, [cosmo]]) in Msun

Give the comoving radius R (in Mpc) in terms of the  Lagrangian mass (in Msun)

# Arguments
- `m`: mass in Msun
- `T`: type of Window (default is `TopHatWindow`)
- `cosmo` : cosmology type (default is `planck18_bkg`)
"""
radius_from_mass(m::Real, ::Type{T} = TopHatWindow, cosmo::BkgCosmology = planck18_bkg) where {T<:Window} = (m / volume_factor(T) / mean_ρ(Matter, 0, cosmo))^(1/3)


dradius_dmass(r::Real, ::Type{T} = TopHatWindow, cosmo::BkgCosmology = planck18_bkg) where {T<:Window} =  r^(-2) / 3.0 / volume_factor(T) / mean_ρ(Matter, 0, cosmo)

# ---------------------------------------------------------
# SMOOTHED VARIANCE OF THE POWER SPECTRUM (AND DERIVATIVES)


@doc raw""" 
   σ2_mps(r, T, matter_ps; kws...)

Give the variance of the function `matter_ps` smoothed on region of size `r` in Mpc.

``{\sigma}^{2}_{\rm mps}(R) = \int_0^{\infty} \mathcal{P}(k) |W(kR)|^2 {\rm d} \ln k``

# Arguments
- `r`: radius in Mpc
- `T`: type of Window 
- `matter_ps`: function (dimensionless power spectrum)
- `kws`: arguments passed to the integration routine `QuadGK.quadgk`
"""
σ2_mps(r::Real, ::Type{T}, matter_ps::Function; kws...) where {T<:Window} = σ2_mps(r, T, matter_ps)
σ2_mps(r::Real, ::Type{TopHatWindow}, matter_ps::Function; kws...)   = log(20.0) - log(r) > -8.0 ? QuadGK.quadgk(lnk -> matter_ps(exp(lnk)) * window_function(exp(lnk) * r, TopHatWindow)^2, -8.0, log(20.0) - log(r), rtol=1e-4; kws...)[1] : 0.0
σ2_mps(r::Real, ::Type{SharpKWindow}, matter_ps::Function; kws...)   = - log(r) > -8.0 ? QuadGK.quadgk(lnk -> matter_ps(exp(lnk)), -8.0, -log(r), rtol=1e-8)[1] : 0.0
σ2_mps(r::Real, ::Type{GaussianWindow}, matter_ps::Function; kws...) = log(4.0) - log(r) > -8.0 ? QuadGK.quadgk(lnk -> matter_ps(exp(lnk)) * window_function(exp(lnk) * r, GaussianWindow)^2, -8.0, log(4.0) - log(r), rtol=1e-4)[1] : 0.0

σ2_mps(r::Real, ::Type{T} = TopHatWindow; cosmology::Cosmology = planck18) where {T<:Window} = σ2_mps(r, T, k->matter_power_spectrum(k, 0.0, cosmology = cosmology, dimensionless = true))
σ_mps(r::Real, ::Type{T} = TopHatWindow; cosmology::Cosmology = planck18) where {T<:Window} = sqrt(σ2_mps(r, T, k->matter_power_spectrum(k, 0.0, cosmology = cosmology, dimensionless = true)))
σ_mps(r::Real, ::Type{T}, matter_ps::Function) where {T<:Window} = sqrt(σ2_mps(r, T, matter_ps))

dσ2_mps_dR(r::Real, ::Type{TopHatWindow}, matter_ps::Function) = log(20.0) - log(r) > -8.0 ? QuadGK.quadgk(lnk -> matter_ps(exp(lnk)) * 2 * window_function(exp(lnk) * r, TopHatWindow) * dwindow_function_dkR(exp(lnk) * r, TopHatWindow) * exp(lnk), -8.0, log(20.0) - log(r), rtol=1e-4)[1] : 0.0
dσ2_mps_dR(r::Real, ::Type{SharpKWindow}, matter_ps::Function) = - matter_ps(1.0 / r) / r
dσ2_mps_dR(r::Real, ::Type{GaussianWindow}, matter_ps::Function) = log(4.0) - log(r) > -8.0 ? QuadGK.quadgk(lnk -> matter_ps(exp(lnk)) * 2 * window_function(exp(lnk) * r, GaussianWindow) * dwindow_function_dkR(exp(lnk) * r, GaussianWindow) * exp(lnk), -8.0, log(4.0) - log(r), rtol=1e-4)[1] : 0.0
dσ2_mps_dR(r::Real, ::Type{T} = TopHatWindow; cosmology::Cosmology = planck18) where {T<:Window} = dσ2_mps_dR(r, T, k->matter_power_spectrum(k, 0.0, cosmology = cosmology, dimensionless = true))
dσ2_mps_dR(r::Real, ::Type{T}, matter_ps::Function) where {T<:Window} = dσ2_mps_dR(r, T, matter_ps)

dσ_mps_dR(r::Real, ::Type{T} = TopHatWindow; cosmology::Cosmology = planck18) where {T<:Window} = dσ_mps_dR(r, T, k->matter_power_spectrum(k, 0.0, cosmology = cosmology, dimensionless = true))
dσ_mps_dR(r::Real, ::Type{T}, matter_ps::Function) where {T<:Window} = 0.5 * dσ2_mps_dR(r, T, matter_ps) / σ_mps(r, T, matter_ps) 

σ2_mps_M(m::Real, ::Type{T}, matter_ps::Function, bkg_cosmology::BkgCosmology = planck18_bkg) where {T<:Window} = σ2_mps(radius_from_mass(m, T, bkg_cosmology), T, matter_ps)
σ_mps_M(m::Real, ::Type{T}, matter_ps::Function, bkg_cosmology::BkgCosmology = planck18_bkg) where {T<:Window} = σ_mps(radius_from_mass(m, T, bkg_cosmology), T, matter_ps)
σ2_mps_M(m::Real, ::Type{T} = TopHatWindow; cosmology::Cosmology = planck18) where {T<:Window} =  σ2_mps(radius_from_mass(m, T, cosmology.bkg), T, cosmology = cosmology)
σ_mps_M(m::Real, ::Type{T} = TopHatWindow; cosmology::Cosmology = planck18) where {T<:Window} =  σ_mps(radius_from_mass(m, T, cosmology.bkg), T, cosmology = cosmology)

function dσ2_mps_dM(m::Real, ::Type{T} = TopHatWindow; cosmology::Cosmology = planck18) where {T<:Window} 
    _r = radius_from_mass(m, T, cosmology.bkg)
    return dσ2_mps_dR(_r, T, cosmology = cosmology) * dradius_dmass(_r, T, cosmology.bkg)
end

function dσ_mps_dM(m::Real, ::Type{T} = TopHatWindow; cosmology::Cosmology = planck18) where {T<:Window} 
    _r = radius_from_mass(m, T, cosmology.bkg)
    return dσ_mps_dR(_r, T, cosmology = cosmology) * dradius_dmass(_r, T, cosmology.bkg)
end

function dσ_dM(m::Real, ::Type{T}, matter_ps::Function, bkg_cosmology::BkgCosmology = planck18_bkg) where {T<:Window} 
    _r = radius_from_mass(m, T, bkg_cosmology)
    return dσ_mps_dR(_r, T, matter_ps) * dradius_dmass(_r, T, bkg_cosmology)
end

function dσ2_mps_dM(m::Real, ::Type{T}, matter_ps::Function, bkg_cosmology::BkgCosmology = planck18_bkg) where {T<:Window} 
    _r = radius_from_mass(m, T, bkg_cosmology)
    return dσ2_mps_dR(_r, T, matter_ps) * dradius_dmass(_r, T, bkg_cosmology)
end
