
const Region5_tup = (
    no = [-0.131_799_836_742_01E2,
    0.685_408_416_344_34E1,
   -0.248_051_489_334_66E-1,
    0.369_015_349_803_33,
   -0.311_613_182_139_25E1,
   -0.329_616_265_389_17],

    Jo = [0,
    1,
    -3,
    -2,
    -1,
    2],

    nr = [0.157_364_048_552_59E-2,
    0.901_537_616_739_44E-3,
    -0.502_700_776_776_48E-2,
    0.224_400_374_094_85E-5,
    -0.411_632_754_534_71E-5,
    0.379_194_548_229_55E-7],

    Ir = [1,
    1,
    1,
    2,
    2,
    3],

    Jr = [1,
    2,
    3,
    3,
    9,
    7],
)
function Region5(Output::Symbol, P, T)
    Pstar = 1.0     #MPa
    Tstar = 1000.0  #K
    no = Region5_tup.no
    Jo = Region5_tup.Jo
    nr = Region5_tup.nr
    Ir = Region5_tup.Ir
    Jr = Region5_tup.Jr



    π = P / Pstar
    τ = Tstar / T

    γo    =  log(π) + sum([no[i]*(τ^Jo[i]) for i=1:6])
    γo_π  =  1/π
    γo_ππ = -1/(π^2)
    γo_τ  =          sum([no[i]*Jo[i]*(τ^(Jo[i]-1)) for i=1:6])
    γo_ττ =          sum([no[i]*Jo[i]*(Jo[i]-1)*(τ^(Jo[i]-2)) for i=1:6])
    γo_πτ = 0

    γr    = sum([nr[i]*(π^Ir[i])*(τ^Jr[i]) for i=1:6])
    γr_π  = sum([nr[i]*Ir[i]*(π^(Ir[i]-1))*(τ^Jr[i]) for i=1:6])
    γr_ππ = sum([nr[i]*Ir[i]*(Ir[i]-1)*(π^(Ir[i]-2))*(τ^Jr[i]) for i=1:6])
    γr_τ  = sum([nr[i]*(π^Ir[i])*Jr[i]*(τ^(Jr[i]-1)) for i=1:6])
    γr_ττ = sum([nr[i]*(π^Ir[i])*Jr[i]*(Jr[i]-1)*(τ^(Jr[i]-2)) for i=1:6])
    γr_πτ = sum([nr[i]*Ir[i]*(π^(Ir[i]-1))*Jr[i]*(τ^(Jr[i]-1)) for i=1:6])

    if Output == :SpecificG            #kJ/kg
        return IF97_R*T*(γo + γr)
    elseif Output == :SpecificF            #kJ/kg
        return IF97_R*T*((γo + γr) - π*(γo_π + γr_π)/1000)
    elseif Output == :SpecificV        #m3/kg
        return IF97_R*T*π*(γo_π + γr_π)/P/1000
    elseif Output == :SpecificU        #kJ/kg
        return IF97_R*T*(τ*(γo_τ + γr_τ) - π*(γo_π + γr_π))
    elseif Output == :SpecificS        #kJ/kgK
        return IF97_R*(τ*(γo_τ + γr_τ) - (γo + γr))
    elseif Output == :SpecificH        #kJ/kg
        return IF97_R*T*τ*(γo_τ + γr_τ)
    elseif Output == :SpecificCP       #kJ/kgK
        return -IF97_R*τ^2*(γo_ττ + γr_ττ)
    elseif Output == :SpecificCV       #kJ/kgK
        return IF97_R*((-τ^2)*(γo_ττ + γr_ττ) - (1 + π*γr_π - τ*π*γr_πτ)^2/(1-(π^2)*γr_ππ))
    elseif Output == :SpeedOfSound     #m/s
        return √(1000*IF97_R*T*(1 + 2*π*γr_π + (π^2)*(γr_π^2))/((1 - (π^2)*(γr_ππ))+((1+π*γr_π-τ*π*γr_πτ)^2)/(τ^2*(γo_ττ + γr_ττ))))
    else
        throw(DomainError(Output, "Unknown value requested."))
    end
end


function mass_gibbs_impl(mt::SinglePT,model::IF97Region{:r5},p,t)
    P = convert_unit(u"Pa",u"MPa",p)
    T = normalize_units(t)
    Pstar = 1.0    #MPa
    Tstar = 1000.0  #K
    no = Region5_tup.no
    Jo = Region5_tup.Jo
    nr = Region5_tup.nr
    Ir = Region5_tup.Ir
    Jr = Region5_tup.Jr

    π = P / Pstar
    τ = Tstar / T
    γo = log(π) + sum(no[i]*τ^Jo[i] for i=1:9)
    γr = sum(nr[i]*π^Ir[i]*(τ-0.5)^Jr[i] for i = 1:6)
    res = IF97_R*T*(γo + γr)
    return convert_unit(u"kJ/kg",u"J/kg",res)
end

function mass_helmholtz_impl(mt::SinglePT,model::IF97Region{:r5},p,t)
    P = convert_unit(u"Pa",u"MPa",p)
    T = normalize_units(t)
    Pstar = 1.0    #MPa
    Tstar = 1000.0  #K
    no = Region5_tup.no
    Jo = Region5_tup.Jo
    nr = Region5_tup.nr
    Ir = Region5_tup.Ir
    Jr = Region5_tup.Jr
    π = P / Pstar
    τ = Tstar / T
    γo    = log(π) + sum(no[i]*τ^Jo[i] for i=1:9)
    γr    = sum(nr[i]*π^Ir[i]*(τ-0.5)^Jr[i] for i = 1:6)
    γo_π  = 1/π
    γr_π  = sum(nr[i]*Ir[i]*π^(Ir[i]-1)*(τ-0.5)^Jr[i] for i=1:6)
    res =  IF97_R*T*(γo + γr - π*(γo_π + γr_π)/1000)
    return convert_unit(u"kJ/kg",u"J/kg",res)
end

function mass_volume_impl(mt::SinglePT,model::IF97Region{:r5},p,t)
    P = convert_unit(u"Pa",u"MPa",p)
    T = normalize_units(t)
    Pstar = 1.0    #MPa
    Tstar = 1000.0  #K
    no = Region5_tup.no
    Jo = Region5_tup.Jo
    nr = Region5_tup.nr
    Ir = Region5_tup.Ir
    Jr = Region5_tup.Jr

    π = P / Pstar
    τ = Tstar / T

    γo_π  = 1/π
    γr_π  = sum(nr[i]*Ir[i]*π^(Ir[i]-1)*(τ-0.5)^Jr[i] for i=1:6)

    res = IF97_R*T*π*(γo_π + γr_π)/P/1000
    return res
end

function mass_internal_energy_impl(mt::SinglePT,model::IF97Region{:r5},p,t)
    P = convert_unit(u"Pa",u"MPa",p)
    T = normalize_units(t)
    Pstar = 1.0    #MPa
    Tstar = 1000.0  #K
    no = Region5_tup.no
    Jo = Region5_tup.Jo
    nr = Region5_tup.nr
    Ir = Region5_tup.Ir
    Jr = Region5_tup.Jr

    π = P / Pstar
    τ = Tstar / T
    γo_τ  = sum(no[i]*Jo[i]*τ^(Jo[i]-1) for i=1:9)
    γr_τ  = sum(nr[i]*π^Ir[i]*Jr[i]*(τ-0.5)^(Jr[i]-1) for i=1:6)
    γo_π  = 1/π
    γr_π  = sum(nr[i]*Ir[i]*π^(Ir[i]-1)*(τ-0.5)^Jr[i] for i=1:6)

    res = IF97_R*T*(τ*(γo_τ + γr_τ) - π*(γo_π + γr_π))

    return convert_unit(u"kJ/kg",u"J/kg",res)
end

function mass_entropy_impl(mt::SinglePT,model::IF97Region{:r5},p,t)
    P = convert_unit(u"Pa",u"MPa",p)
    T = normalize_units(t)
    Pstar = 1.0    #MPa
    Tstar = 1000.0  #K
    no = Region5_tup.no
    Jo = Region5_tup.Jo
    nr = Region5_tup.nr
    Ir = Region5_tup.Ir
    Jr = Region5_tup.Jr

    π = P / Pstar
    τ = Tstar / T
    γo_τ  = sum(no[i]*Jo[i]*τ^(Jo[i]-1) for i=1:9)
    γr_τ  = sum(nr[i]*π^Ir[i]*Jr[i]*(τ-0.5)^(Jr[i]-1) for i=1:6)

    γo    = log(π) + sum(no[i]*τ^Jo[i] for i=1:9)
    γr    = sum(nr[i]*π^Ir[i]*(τ-0.5)^Jr[i] for i = 1:6)

    res = IF97_R*(τ*(γo_τ + γr_τ) - (γo + γr))

    return convert_unit(u"kJ/(kg*K)",u"J/(kg*K)",res)
end

function mass_enthalpy_impl(mt::SinglePT,model::IF97Region{:r5},p,t)
    P = convert_unit(u"Pa",u"MPa",p)
    T = normalize_units(t)
    Pstar = 1.0    #MPa
    Tstar = 1000.0  #K
    no = Region5_tup.no
    Jo = Region5_tup.Jo
    nr = Region5_tup.nr
    Ir = Region5_tup.Ir
    Jr = Region5_tup.Jr

    π = P / Pstar
    τ = Tstar / T
    γo_τ  = sum(no[i]*Jo[i]*τ^(Jo[i]-1) for i=1:9)
    γr_τ  = sum(nr[i]*π^Ir[i]*Jr[i]*(τ-0.5)^(Jr[i]-1) for i=1:6)

    res = IF97_R*T*τ*(γo_τ + γr_τ)
    return convert_unit(u"kJ/kg",u"J/kg",res)
end

function mass_cp_impl(mt::SinglePT,model::IF97Region{:r5},p,t)
    P = convert_unit(u"Pa",u"MPa",p)
    T = normalize_units(t)
       Pstar = 1.0    #MPa
    Tstar = 1000.0  #K
    no = Region5_tup.no
    Jo = Region5_tup.Jo
    nr = Region5_tup.nr
    Ir = Region5_tup.Ir
    Jr = Region5_tup.Jr

    π = P / Pstar
    τ = Tstar / T
    γo_ττ = sum(no[i]*Jo[i]*(Jo[i]-1)*τ^(Jo[i]-2) for i=1:9)
    γr_ττ = sum(nr[i]*π^Ir[i]*Jr[i]*(Jr[i]-1)*(τ-0.5)^(Jr[i]-2) for i=1:6)

    res = -IF97_R*τ^2*(γo_ττ + γr_ττ)
   
    return convert_unit(u"kJ/(kg*K)",u"J/(kg*K)",res)

end

function mass_cv_impl(mt::SinglePT,model::IF97Region{:r5},p,t)
    P = convert_unit(u"Pa",u"MPa",p)
    T = normalize_units(t)
    Pstar = 1.0    #MPa
    Tstar = 1000.0  #K
    no = Region5_tup.no
    Jo = Region5_tup.Jo
    nr = Region5_tup.nr
    Ir = Region5_tup.Ir
    Jr = Region5_tup.Jr

    π = P / Pstar
    τ = Tstar / T
    
    γo_ττ = sum(no[i]*Jo[i]*(Jo[i]-1)*τ^(Jo[i]-2) for i=1:9)
    γr_ττ = sum(nr[i]*π^Ir[i]*Jr[i]*(Jr[i]-1)*(τ-0.5)^(Jr[i]-2) for i=1:6)
    γr_π  = sum(nr[i]*Ir[i]*π^(Ir[i]-1)*(τ-0.5)^Jr[i] for i=1:6)
    γr_πτ = sum(nr[i]*Ir[i]*π^(Ir[i]-1)*Jr[i]*(τ-0.5)^(Jr[i]-1) for i = 1:6)
    γr_ππ = sum(nr[i]*Ir[i]*(Ir[i]-1)*π^(Ir[i]-2)*(τ-0.5)^Jr[i] for i=1:6) 

    res = IF97_R*(-τ^2*(γo_ττ+γr_ττ) - (1+π*γr_π-τ*π*γr_πτ)^2/(1-π^2*γr_ππ))

    return convert_unit(u"kJ/(kg*K)",u"J/(kg*K)",res)
end

function sound_speed_impl(mt::SinglePT,model::IF97Region{:r5},p,t)
    P = convert_unit(u"Pa",u"MPa",p)
    T = normalize_units(t)
    Pstar = 1.0    #MPa
    Tstar = 1000.0  #K
    no = Region5_tup.no
    Jo = Region5_tup.Jo
    nr = Region5_tup.nr
    Ir = Region5_tup.Ir
    Jr = Region5_tup.Jr

    π = P / Pstar
    τ = Tstar / T


    γr_π  = sum(nr[i]*Ir[i]*π^(Ir[i]-1)*(τ-0.5)^Jr[i] for i=1:6)
    γr_ππ = sum(nr[i]*Ir[i]*(Ir[i]-1)*π^(Ir[i]-2)*(τ-0.5)^Jr[i] for i=1:6)
    γr_πτ = sum(nr[i]*Ir[i]*π^(Ir[i]-1)*Jr[i]*(τ-0.5)^(Jr[i]-1) for i = 1:6)
    γo_ττ = sum(no[i]*Jo[i]*(Jo[i]-1)*τ^(Jo[i]-2) for i=1:9)
    γr_ττ = sum(nr[i]*π^Ir[i]*Jr[i]*(Jr[i]-1)*(τ-0.5)^(Jr[i]-2) for i=1:6)

    res =  sqrt((1000*IF97_R*T)*(1+2*π*γr_π+π^2*γr_π^2)/
    ((1-π^2*γr_ππ)+(1+π*γr_π-τ*π*γr_πτ)^2/τ^2/(γo_ττ + γr_ττ)))
    return res
end

