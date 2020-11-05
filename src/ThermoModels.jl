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
using Unitful
using Unitful: @u_str

import Roots
using Roots: find_zeros
using Roots: find_zero

import NLSolversBase
using NLSolversBase: only_fg!,OnceDifferentiable

import LinearAlgebra
using LinearAlgebra: dot,norm

import BangBang
using BangBang: @!

import Optim
using Optim: optimize
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

include("models/cubic.jl")


include("models/iapws95.jl")
include("models/gerg2008.jl")
include("models/lennon2000air.jl")



export IAPWS95, GERG2008
export Saturation 
export equilibria,kvalues
export fugacity_coeff_impl,kvalues,flash_impl,normalizefrac!!

end
