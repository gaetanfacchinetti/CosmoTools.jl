########################################################
# Subfile of the module Cosmojuly.jl
#
# Contains functions related to the transfer functions
#
# author: Gaetan Facchinetti
# email: gaetanfacc@gmail.com
#
########################################################



export transfer_function, TrivialTF, EH98, TransferFunctionModel, EH98_planck18

abstract type TransferFunctionModel end

struct EH98{T<:Real} <: TransferFunctionModel
    
    # Parameters directly related to the BkgCosmology
    Ω_m0_h2::T
    Ω_b0_h2::T
    Ω_χ0_h2::T
    Θ27::T

    # Derived parameters
    z_drag::T
    sound_horizon_Mpc::T
    α_c::T
    α_b::T
    β_c::T
    β_b::T
    k_Silk::T

end

g_func(y::Real) = (-6.0*sqrt(1.0+y) + (2.0+3.0*y)*log((sqrt(1.0+y)+1)/(sqrt(1.0+y)-1.0))) * y


function EH98(Ω_m0_h2::Real, Ω_b0_h2::Real, Ω_χ0_h2::Real, z_eq_mr::Real, k_eq_mr::Real, ::Type{T}; T0_CMB_K::Real = 2.72548) where {T<:Real}
    
    Θ27::T = T0_CMB_K / 2.7

    # z_drag
    b1 = 0.313 * Ω_m0_h2^(-0.419) * (1. + 0.607 * Ω_m0_h2^0.674)
    b2 = 0.238 * Ω_m0_h2^0.223
    z_drag::T = 1291. * Ω_m0_h2^0.251 / (1. + 0.659 * Ω_m0_h2^0.828)  * (1. + b1 * Ω_b0_h2^b2)
    
    # sound horizon
    R_drag = 31.5 * Ω_b0_h2 * Θ27^(-4) * 1e+3 / z_drag
    R_eq   = 31.5 * Ω_b0_h2 * Θ27^(-4) * 1e+3 / z_eq_mr
    sound_horizon_Mpc::T = 2. / (3. * k_eq_mr) * sqrt(6. / R_eq) * log((sqrt(1. + R_drag) + sqrt(R_drag + R_eq))/(1+sqrt(R_eq)))

    # α_c
    a1  = (46.9 * Ω_m0_h2)^0.670 * (1.0 + (32.1 * Ω_m0_h2)^(-0.532) )
    a2  = (12.0 * Ω_m0_h2)^0.424 * (1.0 + (45.0 * Ω_m0_h2)^(-0.582) )
    α_c::T =  a1^(-Ω_b0_h2/Ω_m0_h2) * a2^(-(Ω_b0_h2/Ω_m0_h2)^3) 
    
    # β_c
    b1_2 = 0.944 / (1+(458.0*Ω_m0_h2)^(-0.708))
    b2_2 = (0.395 * Ω_m0_h2)^(-0.0266)
    β_c::T  = 1.0 / (1.0 + b1_2 * ((Ω_χ0_h2/Ω_m0_h2)^b2_2 - 1.0))

    # k_silk, α_b, and β_b
    k_Silk::T = 1.6 * ( Ω_b0_h2^0.52 ) * (Ω_m0_h2^0.73 ) * ( 1.0 + (10.4 * Ω_m0_h2)^(-0.95) )
    α_b::T = 2.07 * k_eq_mr * sound_horizon_Mpc * (1.0 + R_drag)^(-3.0/4.0) * g_func( (1.0 + z_eq_mr) / (1.0 + z_drag) )
    β_b::T = 0.5 + Ω_b0_h2 / Ω_m0_h2 + (3.0 - 2.0 * Ω_b0_h2/Ω_m0_h2) * sqrt( 1.0 + (17.2 * Ω_m0_h2)^2 )

    EH98(convert(T, Ω_m0_h2), convert(T, Ω_b0_h2), convert(T, Ω_χ0_h2), convert(T, Θ27), z_drag, sound_horizon_Mpc, α_c, α_b, β_c, β_b, k_Silk)
end


function EH98(bkg_cosmo::FLRW{<:Real}, ::Type{T} = Float64) where {T<:Real}

    Ω_m0_h2::T = bkg_cosmo.Ω_m0 * bkg_cosmo.h^2
    Ω_b0_h2::T = bkg_cosmo.Ω_b0 * bkg_cosmo.h^2
    Ω_χ0_h2::T = bkg_cosmo.Ω_χ0 * bkg_cosmo.h^2

    z_eq = z_eq_mr(bkg_cosmo)
    k_eq_Mpc = k_eq_mr(bkg_cosmo)

    return EH98(Ω_m0_h2, Ω_b0_h2, Ω_χ0_h2, z_eq, k_eq_Mpc, T, T0_CMB_K = bkg_cosmo.T0_CMB_K)
end

const EH98_planck18 = EH98(planck18_bkg)

function transfer_0_tilde(q::Real, α_c::Real, β_c::Real)::Real
    C = 14.2 / α_c + 386.0 / (1.0 + 69.9 * q^1.08 )
    return log(exp(1.0) + 1.8 * β_c * q ) / (log(exp(1.0) + 1.8 * β_c * q) + C * q^2)
end 

shape_parameter(k::Real, T0_CMB_K::Real, Ω_m0_h2::Real)=  k * (T0_CMB_K/2.7)^2 / Ω_m0_h2
shape_parameter(k::Real, p::EH98)::Real = shape_parameter(k, 2.7 * p.Θ27, p.Ω_m0_h2)

function transfer_cdm(k::Real, p::EH98)::Real
    q = shape_parameter(k, p)
    f = 1.0 / (1.0  +  (k * p.sound_horizon_Mpc / 5.4)^4)
    return f * transfer_0_tilde(q, 1.0, p.β_c) + (1.0 - f) * transfer_0_tilde(q, p.α_c, p.β_c)
end

function s_tilde_Mpc(k::Real, p::EH98)::Real
    β_node = 8.41 * (p.Ω_m0_h2)^0.435
    return p.sound_horizon_Mpc * (1.0 + (β_node/(k * p.sound_horizon_Mpc))^3 )^(-1.0/3.0)
end

function transfer_baryons(k::Real, p::EH98)::Real
    q = shape_parameter(k, p)
    j0 = sinc(k * s_tilde_Mpc(k, p) / π)
    ks = k * p.sound_horizon_Mpc
    return j0 * (transfer_0_tilde(q, 1.0, 1.0) / (1.0 + (ks /5.2)^2 ) + p.α_b * exp(-(k/p.k_Silk)^1.4) / (1.0 + (p.β_b/ks)^3 ) )
end



@doc raw"""
    transfer_function(k, [bkg_cosmo::BkgCosmology], with_baryons = true)

If bkg_cosmo is given computes the transfer function for

"""
function transfer_function(k::Real, bkg_cosmo::FLRW{<:Real}; with_baryons::Bool = true)
    p = EH98(bkg_cosmo)
    tc = transfer_cdm(k, p)
    tb = transfer_baryons(k, p)
    return with_baryons ? p.Ω_b0_h2 / p.Ω_m0_h2 * tb + p.Ω_χ0_h2/p.Ω_m0_h2 * tc : p.Ω_χ0_h2 / p.Ω_m0_h2 * tc
end

function transfer_function(k::Real, p::EH98 = EH98_planck18; with_baryons::Bool = true)
    tc = transfer_cdm(k, p)
    tb = transfer_baryons(k, p)
    return with_baryons ? p.Ω_b0_h2 / p.Ω_m0_h2 * tb + p.Ω_χ0_h2/p.Ω_m0_h2 * tc : p.Ω_χ0_h2 / p.Ω_m0_h2 * tc
end

## Define here the trivial transfer function
struct TrivialTF <: TransferFunctionModel end
TrivialTF(bkg_cosmo::BkgCosmology) = TrivialTF()
transfer_function(k::Real, p::TrivialTF; with_baryons::Bool = true) = 1

## Other transfer functions can be implemented down here

