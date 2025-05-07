# CosmoTools

[![Build Status](https://github.com/gaetanfacchinetti/CosmoTools.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/gaetanfacchinetti/CosmoTools.jl/actions/workflows/CI.yml?query=branch%3Amain)

This is a lightweight Julia module to compute basic cosmological quantities.

## Installation guide

### Recomanded installation

The recomanded way to install CosmoTools is to use the _local_ registry **CosmoRegistry** following the few steps below in a julia shell
```julia
using Pkg

# add the registry to the your own
Pkg.Registry.add(url = "https://github.com/gaetanfacchinetti/CosmoRegistry.git")

# fetch CosmoTools in the local registry
Pkg.add("CosmoTools")

# import CosmoTools
using CosmoTools
```

### For developpers

To develop CosmoTools, you can clone this github repository

```bash
git clone https://github.com/gaetanfacchinetti/CosmoTools.git
```
and install it in a julia shell
```julia
using Pkg
Pkg.develop(path = "<path>/CosmoTools.jl")

using CosmoTools
```

## User manual

The fundamental brick of the CosmoTools is the **BkgCosmology**. You can for instance define a **FLRW** child as shown below.
```julia
using CosmoTools

# Constructor of FLRW cosmology is
# FLRW(h::Real, Ω_χ0::Real, Ω_b0::Real, Ω_k0::Real=0; T0_CMB_K::Real = 2.72548, Neff::Real = 3.04, EdS::Bool = false)
cosmo_bkg::FLRW = FLRW(0.6736, 3.0, 0.0)

# predefined cosmology exists such as
const planck18_bkg::FLRW = FLRW(0.6736, 0.26447, 0.04930)
const edsPlanck18_bkg::FLRW  = FLRW(0.6736, 0.3, 0, EdS = true)
```
A multitude of functions are then available to compute quantities related to this background cosmology.
```julia
# critical density of the Universe
ρ_critical(0.0, cosmo_bkg) 

# age of the universe
universe_age(0.0, cosmo_bkg)
```
In addition, CosmoTools can be used to evaluate the matter power spectrum with a predefined transfer function and even evaluate the halo mass function from the Press-Schechter formalism.
```julia
# define a cosmology from a background cosmology, a primodial curvature power spectrum and a transfer function
# a predifined cosmology is, for instance, that from Planck18 
const planck18::Cosmology = Cosmology("PLANCK18", planck18_bkg, k->power_spectrum_ΛCDM(k, 1e-10*exp(3.044), 0.9649), EH98_planck18)

# the Press-Schechter halo mass function for a top-hat window function in planck18 cosmology at z = 0 is obtained from 
m = 10.0.^range(-5, 16, 200)
dndm = dn_dM.(m, 0, TopHatWindow, PressSchechter; cosmology = planck18, growth_function=growth_factor_Carroll, δ_c = 1.686)
```
## More information
Any comment and contribution is welcome!