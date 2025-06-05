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
# Contains functions related to the transfer functions
#
# author: Gaetan Facchinetti
# email: gaetanfacc@gmail.com
#
##################################################################################


export G_NEWTON, C_LIGHT, MPC_TO_KM, KM_TO_MPC, YR_TO_S, MYR_TO_S
export Unit, Length, Mass, Time, MegaParsecs, KiloMeters, Meters, Msun, MegaYears, Years, Seconds
export constant_G_NEWTON, constant_C_LIGHT, convert_lengths, convert_times

abstract type Unit end
abstract type Length <: Unit end
abstract type Mass <: Unit end
abstract type Time <: Unit end

abstract type MegaParsecs <: Length end
abstract type KiloMeters <: Length end
abstract type Meters <: Length end

abstract type Msun <: Mass end

abstract type MegaYears <: Time end
abstract type Years <: Time end
abstract type Seconds <: Time end

# First define some basic constants (the only one we will need)
const G_NEWTON::Float64 = 4.517103049894965e-48 # in Mpc^3 / Msun / s^(2)
const C_LIGHT::Float64 = 9.715611890180198e-15 # Mpc / s 
const MPC_TO_KM::Float64 = 3.085677581491367e19
const KM_TO_MPC::Float64 = 1.0/MPC_TO_KM
const YR_TO_S::Float64 = 31557600 
const MYR_TO_S::Float64 = 31557600e+6


constant_G_NEWTON(::Type{KiloMeters},  ::Type{Msun}, ::Type{Seconds}, ::Type{Float64})::Float64 = 1.3271244e11
constant_G_NEWTON(::Type{KiloMeters},  ::Type{Msun}, ::Type{Seconds}, ::Type{Float32})::Float32 = 1.3271244f11
constant_G_NEWTON(::Type{MegaParsecs}, ::Type{Msun}, ::Type{Seconds}, ::Type{Float64})::Float64 = 4.517103049894965e-48
constant_G_NEWTON(::Type{MegaParsecs}, ::Type{Msun}, ::Type{Seconds}, ::Type{Float32}) = throw("Impossible to define G_NEWTON in these units with Float32")

constant_G_NEWTON(::Type{L} = MegaParsecs, ::Type{M} = Msun, ::Type{T} = Seconds, ::Type{S} = Float64) where {L<:Length, M<:Mass, T<:Time, S<:AbstractFloat} = constant_G_NEWTON(L, M, T, S)

constant_C_LIGHT(::Type{MegaParsecs}, ::Type{Seconds}, ::Type{Float64})::Float64 = 9.715611890180198e-15
constant_C_LIGHT(::Type{MegaParsecs}, ::Type{Seconds}, ::Type{Float32})::Float32 = 9.715612f-15
constant_C_LIGHT(::Type{KiloMeters}, ::Type{Seconds}, ::Type{Float64})::Float64  = 2.99792458e+5
constant_C_LIGHT(::Type{KiloMeters}, ::Type{Seconds}, ::Type{Float32})::Float32  = 2.99792458f+5
constant_C_LIGHT(::Type{Meters}, ::Type{Seconds}, ::Type{Float64})::Float64  = 2.99792458e+8
constant_C_LIGHT(::Type{Meters}, ::Type{Seconds}, ::Type{Float32})::Float32  = 2.99792458f+8

constant_C_LIGHT(::Type{L} = MegaParsecs, ::Type{T} = Seconds, ::Type{S} = Float64) where {L<:Length, T<:Time, S<:AbstractFloat} = constant_C_LIGHT(L, T, S)

convert_lengths(::Type{MegaParsecs}, ::Type{KiloMeters}, ::Type{Float64})::Float64 = 3.085677581491367e19
convert_lengths(::Type{MegaParsecs}, ::Type{KiloMeters}, ::Type{Float32})::Float32 = 3.0856776f19
convert_lengths(::Type{KiloMeters}, ::Type{MegaParsecs}, ::Type{Float64})::Float64 = 3.2407792894443654e-20
convert_lengths(::Type{KiloMeters}, ::Type{MegaParsecs}, ::Type{Float32})::Float32 = 3.2407793f-20

convert_lengths(::Type{L1}, ::Type{L2}, ::Type{S} = Float64) where {L1<:Length, L2<:Length, S<:AbstractFloat} = convert_lengths(L1, L2, S)
convert_lengths(length::S, ::Type{L1}, ::Type{L2}) where {L1<:Length, L2<:Length, S<:AbstractFloat} = length * convert_lengths(L1, L2, S)

convert_times(::Type{Years}, ::Type{Seconds}, ::Type{Float64})::Float64 = 31557600.0
convert_times(::Type{Years}, ::Type{Seconds}, ::Type{Float32})::Float32 = 31557600.0f0
convert_times(::Type{MegaYears}, ::Type{Seconds}, ::Type{Float64})::Float64 = 31557600e6
convert_times(::Type{MegaYears}, ::Type{Seconds}, ::Type{Float32})::Float32 = 31557600.0f6

convert_times(::Type{T1}, ::Type{T2}, ::Type{S} = Float64) where {T1<:Time, T2<:Time, S<:AbstractFloat} = convert_times(T1, T2, S)
convert_times(time::S, ::Type{T1}, ::Type{T2}) where {T1<:Time, T2<:Time, S<:AbstractFloat} = time * convert_times(T1, T2, S)