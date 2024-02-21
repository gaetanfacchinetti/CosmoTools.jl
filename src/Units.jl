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

# First define some basic constants (the only one we will need)
const G_NEWTON::Float64 = 4.517103049894965e-48 # in Mpc^3 / Msun / s^(2)
const C_LIGHT::Float64 = 9.715611890180198e-15 # Mpc / s 
const MPC_TO_KM::Float64 = 3.085677581491367e19
const KM_TO_MPC::Float64 = 1.0/MPC_TO_KM
const YR_TO_S::Float64 = 31557600 
const MYR_TO_S::Float64 = 31557600e+6

