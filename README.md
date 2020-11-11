# ThermoModels

This package aims to provide as many thermodynamic models as posible, including equlibria solvers and AD calculation of properties. uses ForwardDiff for derivation, Roots.jl for root-finding, Unitful.jl for units and ThermoState.jl for state specification.
## Installation

```julia-repl
julia> ]
(@v1.5) pkg> add https://github.com/longemen3000/ThermoModels.jl
```

## Equilibria Models
at the moment, the models included capable of performing a phase equilibria are:

- Van der Waals
- Redlich Kwong
- Soave Redlich Kwong
- Peng Robinson
- GERG2008
- IAPWS95

When models share parameters, you can create a model from the parameters of another model:

```julia
using ThermoModels
gerg_watero2 = GERG2008(:H2O,:O20) #this is the order of the compositions
pr_watero2 = PR(gerg_watero2) #

iapws_water = IAPWS95()
srk_water = SRK(iapws_water)
```

## `equilibria`
thermodynamical equilibria can be called from the `equilibria` function. this function takes a feed state and returns one or more resulting states. The `equilibria` function, and the rest of property functions, follow the proposed ThermoState syntax of `property(model,state,unit)`.
`equilibria` doesnt need units as it returns other states, so the syntax is just `equilibria(model,state)`.

```julia
using ThermoModels, ThermoState,Unitful
st = state(t=60.0u"°C",p = 1u"atm",xn=[0.6,0.4])
#0.6 moles of water, 0.4 moles of O2
gerg_watero2 = GERG2008(:H2O,:O20)
equilibria(gerg_watero2,st)
```
the input `ThermodynamicState` is a Pressure-Temperature-feed input point,so a PT flash is done. the result is more `ThermodynamicState`s:

```julia
(state(mol_v = 1.83213e-5[m^3 mol^-1], t = 333.15[K], xn = [1.0, 2.10949e-7], moles = 0.501819[mol], phase = liquid), 
state(mol_v = 0.0273034[m^3 mol^-1], t = 333.15[K], xn = [0.19708, 
0.80292], moles = 0.498181[mol], phase = gas))
```

as models can be more complicated, an speed advantage can be done by doing the calculations of liquid and gas phases on parallel. That is why ThermoModels uses threads to accelerate computations when possible. disabling threading can be done by passing an options named tuple on the creation of the `ThermodynamicState`

```julia
julia> @btime equilibria($gerg_watero2,$st)
  14.443 ms (10415 allocations: 1.14 
MiB)
opts = (;threaded=false)
st2 = state(t=60.0u"°C",
p = 1u"atm",
xn=[0.6,0.4],
options = opts
)

julia> @btime equilibria($gerg_watero2,$st2)
  17.614 ms (10278 allocations: 1.12 
MiB)
```

At the moment, the following equilibria solvers are available, defined on the states that can be passed as input:

- `state(sat=true,p=p0)`: single component, saturation pressure provided
- `state(sat=true,t=t0)`: single component, saturation temperature provided
- `state(p=p0,t=t0,xn=xn0)`: multicomponent flash,pressure and temperature provided.

every model requires 2 properties to participate on a equilibria:
- `fugacity_coeff_impl`
- `pressure_impl`

in pure helmholtz equations, the pressure and fugacity coeffients are obtained by differenciating the molar helmoltz energy (`mol_helmholtz_impl`), but specialized implementations can be provided. for example, Cubics provide an optimized implementation that reuses and reduces the amount of times that their `a` and `b` parameters are calculated.

For total helmholtz equations, the following equations are provided:

-`pressure`
-`mol_helmholtz`,`mass_helmholtz`,`total_helmholtz`
-`mol_gibbs`,`mass_gibbs`,`total_gibbs`
-`mol_enthalpy`,`mass_enthalpy`,`total_enthalpy`
-`mol_internal_energy`,`mass_internal_energy`,`total_internal_energy`
-`mol_entropy`,`mass_entropy`,`total_entropy`
-`mol_cv`,`mass_cv`
-`sound_speed`

with more to come. all functions accept the following states:

- `state(v=v0,t=t0,...)`: Volume-temperature. calls the primal function, so is pretty fast, but the volume is not interchangeable between different models.
- `state(p=p0,t=t0,...)`: Pressure-temperature.
finds the gas and liquid phases, and returns the one with the least amount of gibbs energy

- `state(p=p0,t=t0,phase=:liquid...)`: Pressure-temperature-Liquid phase.
finds the gas and liquid phases, and returns the liquid volume root

- `state(p=p0,t=t0,phase=:gas...)`: Pressure-temperature-Liquid phase.
finds the gas and liquid phases, and returns the gas volume root

## Saturation Models

There are also saturation models, those models have two properties, `temperature` and `pressure`. and only one possible state for each:

- `temperature(model,state(sat=true,p=psat))`
- `pressure(model,state(sat=true,t=tsat))`

the Saturation models provided by this package are:

- `Saturation.Antoine`
- `Saturation.LeeKesler`
- `Saturation.Edalat`
- `Saturation.Sanjari`
- `Saturation.WaterSat` (IF97 equations for water, for the complete IF97 model, check WaterIF97.jl)
- `Saturation.VdWSatApprox`(Padé aproximant for VdW model)
- `Saturation.RKSatApprox`(Padé aproximant for RK model)
- `Saturation.PRSatApprox`(Padé aproximant for PR model)

Helmholtz models use some of those equations to provide an initial guess to their single phase equilibria.


# WARNING!

this package is in *HEAVY* construction, expect a lot of changes


