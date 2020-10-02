# ThermoModels

This package aims to provide as many thermodynamic models as posible, including equlibria solvers and AD calculation of properties. uses ForwardDiff for derivation, Roots.jl for root-finding, Unitful.jl for units and ThermoState.jl for state specification.

one example of what this package provides:
```julia
Water1 = ThermoModels.GERG2008(:H2O) #water component of GERG2008
Water2 = ThermoModels.IAPWS95() #IAPWS helmholtz model
s100 = state(sat=true,t=100u"°C") #saturated. 100 C

```
lets see the GERG model first
```julia-repl
julia> equilibria(Water1,s100)
(3 Constant Property Specifications:
 Molar volume : 1.8803e-5 m^3 mol^-1
 Temperature : 373.15 K
 Phase specification : liquid, 3 Constant Property Specifications:
 Molar volume : 0.0300984 m^3 mol^-1
 Temperature : 373.15 K
 Phase specification : gas)
```

what about the IAPWS95 one?

```julia-repl
julia> equilibria(Water2,s100)
(3 Constant Property Specifications:
 Molar volume : 1.87982e-5 m^3 mol^-1
 Temperature : 373.15 K
 Phase specification : liquid, 3 Constant Property Specifications:
 Molar volume : 0.0301173 m^3 mol^-1
 Temperature : 373.15 K
 Phase specification : gas)
```

equilibrium calculation are expensive, maybe using an aproximation?

```julia-repl
julia> SatAprox = Saturation.LeeKesler(Water1)      
ThermoModels.LeeKesler{Float64}(647.096, 2.2064e7, 0.344861)

julia> @to_units pressure(SatAprox,s100,u"atm")     
0.9005010371916065 atm
```
Horrible aproximation, as expected from water. lets use a model that doesnt require an expensive iteration scheme to calculate saturations, but is more precise:



```julia-repl
julia> Water3 = WaterIF97()

julia> @to_units pressure(Water3,s100,u"atm")       
1.0009176207383188 atm
```

that's much better!, whats the reduced pressure of that value?

```julia-repl
julia> using ThermoState.StatePoints;
p = pressure(Water3,s100) ;      
pc = pressure(Water3,CriticalPoint());
p/pc

0.004596536345237044
```

nice to know!

# WARNING!

this package is in *HEAVY* construction, expect a lot of changes


