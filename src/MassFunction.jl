
export dn_dM, MassFunctionType, PressSchechter, SethTormen


###########################
# PRESS-SCHECHTER FORMALISM
###########################

abstract type MassFunctionType end
abstract type PressSchechter <: MassFunctionType end
abstract type SethTormen <: MassFunctionType end

f_mass_function(ν::Real, ::Type{PressSchechter}) = sqrt(2.0 / π) * exp(-ν^2 / 2.0)

function f_mass_function(ν::Real, ::Type{SethTormen}) 

    _a = 0.707;
    _νp = sqrt(_a) * ν;
    _q = 0.3;
    _A = 0.3222;

    return _νp / ν * _A * (1.0 + _νp^(-2. * _q)) * sqrt(2.0 * _a / π) * exp(- _νp^2 / 2.)
end

function dn_dM(M_Msun::Real, 
                z::Real = 0,
                ::Type{T} = TopHatWindow,
                ::Type{S} = PressSchechter;
                cosmology::Cosmology = planck18,
                growth_function::Function = growth_factor_Carroll,
                δ_c = 1.686) where {T<:Window, S<:MassFunctionType}

    _D_z = growth_function(z, cosmology.bkg) / growth_function(0.0, cosmology.bkg)
    _σ = σ_mps_M(M_Msun, T, cosmology = cosmology) * _D_z
    _dσ_dM = dσ_mps_dM(M_Msun, T,  cosmology = cosmology) * _D_z
    _ν = δ_c / _σ

    return  mean_ρ(Matter, 0.0, cosmology.bkg) / M_Msun * abs(_dσ_dM) / _σ  * _ν  * f_mass_function(_ν, S)

end

""" Mass function for arbitrary matter power_spectrum """
function dn_dM(M_Msun::Real, 
    z::Real, 
    matter_ps::Function,
    ::Type{T} = TopHatWindow, 
    ::Type{S} = PressSchechter;
    bkg_cosmology::BkgCosmology = planck18_bkg,    
    growth_function::Function = growth_factor_Carroll,
    δ_c = 1.686) where {T<:Window, S<:MassFunctionType}
    
    _D_z = growth_function(z, bkg_cosmology) / growth_function(0.0, bkg_cosmology)
    _σ = σ_mps_M(M_Msun, T, matter_ps, bkg_cosmology) * _D_z
    _dσ_dM = dσ_dM(M_Msun, T, matter_ps, bkg_cosmology) * _D_z
    _ν = δ_c / _σ

    return  mean_ρ(Matter, 0.0, bkg_cosmology) / M_Msun * abs(_dσ_dM) / _σ  * _ν  * f_mass_function(_ν, S)
end


