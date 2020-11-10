module ThermoModels

using ThermoState
using ThermoState.Types
using ThermoState.QuickStates
using ThermoState.StatePoints

import ThermoState: pressure,temperature,mass,moles,molar_mass
import ThermoState: mass_volume, mol_volume, total_volume
import ThermoState: mass_density, mol_density
import ThermoState: mass_enthalpy, mol_enthalpy,total_enthalpy
import ThermoState: mass_gibbs, mol_gibbs, total_gibbs
import ThermoState: mass_helmholtz, mol_helmholtz, total_helmholtz
import ThermoState: mass_internal_energy, mol_internal_energy, total_internal_energy
import ThermoState: mass_entropy, mol_entropy, total_entropy
import ThermoState: mass_fraction, mol_fraction
import ThermoState: mass_number, mol_number
import ThermoState: options, phase, quality
import ThermoState: mol_cp,mol_cv,mass_cv,mass_cp
import ThermoState: sound_speed
import ThermoState: molecular_weight


using StaticArrays
using DiffResults
using ForwardDiff

#units
using Unitful
using Unitful: @u_str

#for cubic solvers
import Polynomials

import Roots
using Roots: find_zero

#for pt-flash
import NLSolversBase
using NLSolversBase: only_fg!,OnceDifferentiable

import Optim
using Optim: optimize


import LinearAlgebra
using LinearAlgebra: dot,norm

#to support static arrays correctly
import BangBang
using BangBang: @!


# Write your package code here.
include("base/utils.jl")
include("base/suva.jl")
include("base/types.jl")
include("base/mixing.jl")
include("base/rr_flash.jl")

#this is for using saturation models on estimation of helmholtz models
include("base/saturation.jl")
include("models/saturation.jl")

include("base/helmholtz/helmholtz.jl")

include("base/helmholtz/cubic.jl")
include("base/helmholtz/volume_solver.jl")

include("base/helmholtz/equilibria_single.jl")
include("base/helmholtz/equilibria_pt.jl")


#each model file exports its own model
include("models/iapws95.jl")
include("models/gerg2008.jl")
include("models/cubic/cubics.jl")

#include("models/lennon2000air.jl")


export VanDerWaals, RedlichKwong
export IAPWS95, GERG2008
export Saturation 
export equilibria,kvalues
export compressibility_factor, acentric_factor
export fugacity_coeff_impl,kvalues,flash_impl,normalizefrac!!

end
