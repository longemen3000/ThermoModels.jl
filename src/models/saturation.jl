

#this includes methods for aproximation for saturation pressures and 
#temperatures
"""
    Antoine(a,b,c,from,to,base)
    Antoine(a,b,c;from =u"K",to = u"Pa",base=ℯ)

## Antoine Vapor Pressure

Calculates vapor pressure of a chemical using the Antoine equation. Parameters `a`, `b`, and `c` are chemical-dependent. Parameters can be found in numerous sources; however units of the coefficients used vary. Originally proposed by Antoine (1888)

The vapor pressure is given by:

    log(base,p) = a + b/(c+t)

by default, `t` in in kelvin, `p` in Pa and the base is natural (ℯ)

the coefficients are normally defined for a specific combination of pressure and temperature units, you can specify those using the `from` (temperature) and `to` (pressure) arguments, passing a corresponding Unitful unit.

## References

1. Poling, Bruce E. The Properties of Gases and Liquids. 5th edition. New York: McGraw-Hill Professional, 2000.
2. Antoine, C. 1888. Tensions des Vapeurs: Nouvelle Relation Entre les Tensions et les Tempé. Compt.Rend. 107:681-684.
3. Yaws, Carl L. The Yaws Handbook of Vapor Pressure: Antoine Coefficients. 1 edition. Houston, Tex: Gulf Publishing Company, 2007.
"""
struct Antoine{V,B,F,T} <: SaturationModel
    a::V
    b::V
    c::V
    base::B
    from::F
    to::T
end

Antoine(a,b,c;from =u"K",to = u"Pa",base=ℯ) = Antoine(a,b,c,base,from,to)

function pressure_impl(mt::SingleSatT,model::T,t) where T<: Antoine
    _t = convert_unit(u"K",model.from,t)
    _a = model.a
    _b = model.b
    _c = model.c
    _base = model.base
    p = _base^(_a+_b/(_c+_t))
    return convert_unit(model.to,u"Pa",p)
end


function temperature_impl(mt::SingleSatP,model::T,p) where T<:Antoine
    _p = convert_unit(u"Pa",model.to,p)
    _a = model.a
    _b = model.b
    _c = model.c
    _base = model.base
    _t = _b/(log(_base,_p)-_a) -_c
    return convert_unit(model.from,u"K",_t)
end
"""
    LeeKesler(;tc,pc,ω)
    LeeKesler(tc,pc,ω)

## Lee Kesler Vapor Pressure

Calculates vapor pressure of a fluid at arbitrary temperatures using a
corresponding states relationship by [1]; requires a chemical's critical temperature and
acentric factor.
    
The vapor pressure is given by:

    ln(Pᵣ) = f⁰ + ωf¹
    f⁰ = 5.92714 - 6.09648/Tᵣ - 1.28862*ln(Tᵣ) + 0.169347*Tᵣ^6
    f¹ = 15.2518 - 15.6875/Tᵣ - 13.4721*ln(Tᵣ) + 0.43577*Tᵣ^6

## References
1. Lee, Byung Ik, and Michael G. Kesler. "A Generalized Thermodynamic Correlation Based on Three-Parameter Corresponding States." AIChE Journal 21, no. 3 (1975): 510-527. doi:10.1002/aic.690210313.
2. Reid, Robert C..; Prausnitz, John M.;; Poling, Bruce E. The Properties of Gases and Liquids. McGraw-Hill Companies, 1987.

"""
struct LeeKesler{T} <: SaturationModel
    tc::T
    pc::T
    ω::T
end



function LeeKesler(;tc,pc,ω)
    _tc = ThermoState.normalize_units(tc)
    _pc = ThermoState.normalize_units(pc)
    return LeeKesler(_tc,_pc,ω)
end

function LeeKesler(model::ThermoModel)
    _tc = temperature(model,CriticalPoint(),u"K")
    _pc = pressure(model,CriticalPoint(),u"Pa")
    _ω = acentric_factor(model)
    return LeeKesler(_tc,_pc,_ω)
end

function temperature(model::LeeKesler,st::CriticalPoint,unit=u"K")
    return convert_unit(u"K",unit,model.tc)
end

function pressure(model::LeeKesler,st::CriticalPoint,unit=u"Pa")
    return convert_unit(u"Pa",unit,model.pc)
end

function acentric_factor(model::LeeKesler)
    return model.ω
end

function pressure_impl(mt::SingleSatT,model::LeeKesler,t)
    Pc = model.pc
    Tᵣ = t/model.tc
    logTᵣ = log(Tᵣ)
    ω = model.ω
    Tᵣ6 = Tᵣ*Tᵣ
    Tᵣ6 *= Tᵣ6*Tᵣ6
    f⁰ = 5.92714 - 6.09648/Tᵣ - 1.28862*logTᵣ+ 0.169347*Tᵣ6
    f¹ = 15.2518 - 15.6875/Tᵣ - 13.4721*logTᵣ + 0.43577*Tᵣ6
    return Pc*exp(muladd(f¹,ω,f⁰))
end
 
initial_temperature(model::LeeKesler,p) = critical_sat_interpolation(model,p)


"""
    Wagner(;tc,pc,a,b,c,d)
    Wagner(tc,pc,a,b,c,d)

# Wagner vapor pressure

Calculates vapor pressure using the Wagner equation (2.5, 5 form).

Requires critical temperature and pressure as well as four coefficients specific to each chemical.

The vapor pressure is given by:

    ln(Pᵣ)=  (aτ  + bτ^1.5 + cτ^2.5+ dτ^5)/Tᵣ

    τ = 1 - Tᵣ

## References

1. Wagner, W. “New Vapour Pressure Measurements for Argon and Nitrogen and a New Method for Establishing Rational Vapour Pressure Equations.” Cryogenics 13, no. 8 (August 1973): 470-82. doi:10.1016/0011-2275(73)90003-9
2. Poling, Bruce E. The Properties of Gases and Liquids. 5th edition. New York: McGraw-Hill Professional, 2000.

"""
struct Wagner{T} <: SaturationModel
    tc::T
    pc::T
    a::T
    b::T 
    c::T 
    d::T 
end

function Wagner(;tc,pc,a,b,c,d)
    _tc = normalize_units(tc)
    _pc = normalize_units(pc)
    return Wagner(_tc,_pc,a,b,c,d)
end

initial_temperature(model::Wagner,p) = critical_sat_interpolation(model,p)

function temperature(model::Wagner,st::CriticalPoint,unit=u"K")
    return convert_unit(u"K",unit,model.tc)
end

function pressure(model::Wagner,st::CriticalPoint,unit=u"Pa")
    return convert_unit(u"Pa",unit,model.pc)
end


function pressure_impl(mt::SingleSatT,model::Wagner,t)
    Tc = model.tc
    Pc = model.pc
    Tr = t/Tc
    τ = one(Tr) - Tr
    τhalf = sqrt(tau)
    bτ = τ*τhalf
    cτ = bτ*τ
    dτ = cτ*cτ
    return Pc*exp((a*τ + b*bτ+ c*cτ + d*dτ)/Tr)
end

"""
WagnerOriginal(;tc,pc,a,b,c,d)
WagnerOriginal(tc,pc,a,b,c,d)

# Wagner (original) vapor pressure

Calculates vapor pressure using the Wagner equation (3, 6 form).

Requires critical temperature and pressure as well as four coefficients specific to each chemical.

The vapor pressure is given by:

    ln(Pᵣ)=  (aτ  + bτ^1.5 + cτ^3 dτ^6)/Tᵣ

    τ = 1 - Tᵣ

## References

1. Poling, Bruce E. The Properties of Gases and Liquids. 5th edition. New York: McGraw-Hill Professional, 2000.
2. McGarry, Jack. “Correlation and Prediction of the Vapor Pressures of Pure Liquids over Large Pressure Ranges.” Industrial & Engineering Chemistry Process Design and Development 22, no. 2 (April 1, 1983): 313-22. doi:10.1021/i200021a023.
"""
struct WagnerOriginal{T} <: SaturationModel
    tc::T
    pc::T
    a::T
    b::T 
    c::T 
    d::T 
end

function WagnerOriginal(;tc,pc,a,b,c,d)
    _tc = normalize_units(tc)
    _pc = normalize_units(pc)
    return WagnerOriginal(_tc,_pc,a,b,c,d)
end

function temperature(model::WagnerOriginal,st::CriticalPoint,unit=u"K")
    return convert_unit(u"K",unit,model.tc)
end

function pressure(model::WagnerOriginal,st::CriticalPoint,unit=u"Pa")
    return convert_unit(u"Pa",unit,model.pc)
end


function pressure_impl(mt::SingleSatT,model::WagnerOriginal,t)
    Tc = model.tc
    Pc = model.pc
    Tr = t/Tc
    τ = one(Tr) - Tr
    τhalf = sqrt(tau)
    bτ = τ*τhalf
    cτ = bτ*bτ
    dτ = cτ*cτ
    return Pc*exp((a*τ + b*bτ+ c*cτ + d*dτ)/Tr)
end

initial_temperature(model::WagnerOriginal,p) = critical_sat_interpolation(model,p)

"""
    TRCAntoine(tc,t0,A,B,C,n,E,F)
    TRCAntoine(;tc,t0,B,A,C,n,E,F)

# TRC extended Antoine model for vapor pressure

Calculates vapor pressure of a chemical using the TRC Extended Antoine equation. Parameters are chemical dependent, and said to be from the Thermodynamics Research Center (TRC) at Texas A&M. Coefficients for various chemicals can be found in [1].

The vapor pressure is given by:

    log10(Psat) = A - B/(t + C) + 0.43429*x^n + Ex^8 + Fx^12

    x = max((t-t₀-273.15)/tc, 0)

## References

1. Poling, Bruce E. The Properties of Gases and Liquids. 5th edition. New York: McGraw-Hill Professional, 2000.
"""

struct TRCAntoine{T} <: SaturationModel
    tc::T
    t0::T
    A::T
    B::T 
    C::T 
    n::T
    E::T
    F::T
end

TRCAntoine(;tc,t0,B,A,C,n,E,F) = TRCAntoine(tc,t0,A,B,C,n,E,F)

function temperature(model::TRCAntoine,st::CriticalPoint,unit=u"K")
    return convert_unit(u"K",unit,model.tc)
end

function pressure_impl(mt::SingleSatT,model::TRCAntoine,t)
    _x = (t - to - 273.15)/Tc
    x = max(_x,zero(_x))
    A = model.A
    B = model.B
    C = model.C
    n = model.n
    E = model.E
    F = model.F
    x4 = x^4
    return exp10(A - B/(t+C) + 0.43429*x^n + x4*x4*(E + F*x4))
end

function initial_temperature(model::TRCAntoine,p)
    A = model.A
    B = model.B
    C = model.C
    return B/(log10(p)-A) -C
end


"""
    AmbroseWalton(;tc,pc,ω)
    AmbroseWalton(tc,pc,ω)

# Ambrose-Walton relation for vapor pressure

Calculates vapor pressure of a fluid at arbitrary temperatures using a corresponding states relationship by [1]; requires a chemical’s critical temperature and acentric factor.

The vapor pressure is given by:

    ln(Pᵣ) = f⁰ + f¹ω + f²ω^2
    f⁰ = -5.97616τ+ 1.29874τ^1.5 - 0.60394τ^2.5 - 1.06841τ^5
    f¹ = -5.03365τ + 1.11505τ^1.5 - 5.41217τ^2.5 - 7.46628τ^5
    f² = -0.64771τ + 2.41539τ^1.5 - 4.26979τ^2.5 + 3.25259τ^5
    τ = 1 - Tᵣ

## References

1. Ambrose, D., and J. Walton. “Vapour Pressures up to Their Critical Temperatures of Normal Alkanes and 1-Alkanols.” Pure and Applied Chemistry 61, no. 8 (1989): 1395-1403. doi:10.1351/pac198961081395.
2. Poling, Bruce E. The Properties of Gases and Liquids. 5th edition. New York: McGraw-Hill Professional, 2000.
"""

struct AmbroseWalton{T} <: SaturationModel
    tc::T
    pc::T
    ω::T
end

function AmbroseWalton(;tc,pc,ω)
    _tc = ThermoState.normalize_units(tc)
    _pc = ThermoState.normalize_units(pc)
    return AmbroseWalton(_tc,_pc,ω)
end

function AmbroseWalton(model::ThermoModel)
    _tc = temperature(model,CriticalPoint(),u"K")
    _pc = pressure(model,CriticalPoint(),u"Pa")
    _ω = acentric_factor(model)
    return AmbroseWalton(_tc,_pc,_ω)
end

function temperature(model::AmbroseWalton,st::CriticalPoint,unit=u"K")
    return convert_unit(u"K",unit,model.tc)
end

function pressure(model::AmbroseWalton,st::CriticalPoint,unit=u"Pa")
    return convert_unit(u"Pa",unit,model.pc)
end

function acentric_factor(model::AmbroseWalton)
    return model.ω
end
function pressure_impl(mt::SingleSatT,model::AmbroseWalton,_t)
    t = normalize_units(_t)
    Pc = model.pc
    τ = 1 - t/model.tc
    ω = model.ω
    τ15 = sqrt(τ)*τ
    τ25 = τ15*τ
    τ5 = τ25*τ25
    f⁰ = -5.97616τ+ 1.29874τ15 - 0.60394τ25 - 1.06841τ5
    f¹ = -5.03365τ + 1.11505τ15 - 5.41217τ25 - 7.46628τ5
    f² = -0.64771τ + 2.41539τ15 - 4.26979τ25 + 3.25259τ5
    
    return Pc*exp(muladd(muladd(f²,ω,f¹),ω,f⁰))
end

initial_temperature(model::AmbroseWalton,p) = critical_sat_interpolation(model,p)

"""
    Sanjari(;tc,pc,ω)
    Sanjari(tc,pc,ω)

# Sanjari relation for vapor pressure

Calculates vapor pressure of a fluid at arbitrary temperatures using a corresponding states relationship by [1]; requires a chemical’s critical temperature and acentric factor. Although developed for refrigerants, this model should have some general predictive ability.

The vapor pressure is given by:

    ln(Pᵣ) = f⁰ + f¹ω + f²ω^2
    f⁰ = 6.83377 + -5.76051/Tᵣ + 0.90654*log(Tᵣ) + -1.16906*Tᵣ^1.9
    f¹ = 5.32034 + -28.1460/Tᵣ + -58.0352*log(Tᵣ) + 23.57466*Tᵣ^1.9
    f² = 18.19967 + 16.33839/Tᵣ + 65.6995*log(Tᵣ) + -35.9739*Tᵣ^1.9


## References

(1, 2) Sanjari, Ehsan, Mehrdad Honarmand, Hamidreza Badihi, and Ali Ghaheri. “An Accurate Generalized Model for Predict Vapor Pressure of Refrigerants.” International Journal of Refrigeration 36, no. 4 (June 2013): 1327-32. doi:10.1016/j.ijrefrig.2013.01.007.
"""

struct Sanjari{T} <: SaturationModel
    tc::T
    pc::T
    ω::T
end

function Sanjari(;tc,pc,ω)
    _tc = ThermoState.normalize_units(tc)
    _pc = ThermoState.normalize_units(pc)
    return Sanjari(_tc,_pc,ω)
end

function Sanjari(model::ThermoModel)
    _tc = temperature(model,CriticalPoint(),u"K")
    _pc = pressure(model,CriticalPoint(),u"Pa")
    _ω = acentric_factor(model)
    return Sanjari(_tc,_pc,_ω)
end

function temperature(model::Sanjari,st::CriticalPoint,unit=u"K")
    return convert_unit(u"K",unit,model.tc)
end

function pressure(model::Sanjari,st::CriticalPoint,unit=u"Pa")
    return convert_unit(u"Pa",unit,model.pc)
end

function acentric_factor(model::Sanjari)
    return model.ω
end

function pressure_impl(mt::SingleSatT,model::Sanjari,_t)
    t = normalize_units(_t)
    ω = model.ω
    Tr = t/model.tc
    Tr_inv = 1.0/Tr
    log_Tr = log(Tr)
    Tr_19 = Tr^1.9
    f⁰ = 6.83377 + -5.76051*Tr_inv + 0.90654*log_Tr + -1.16906*Tr_19
    f¹ = 5.32034 + -28.1460*Tr_inv + -58.0352*log_Tr + 23.57466*Tr_19
    f² = 18.19967 + 16.33839*Tr_inv + 65.6995*log_Tr + -35.9739*Tr_19
    

    return model.pc*exp(muladd(muladd(f²,ω,f¹),ω,f⁰))
end



initial_temperature(model::Sanjari,p) = critical_sat_interpolation(model,p)


"""
    Edalat(;tc,pc,ω)
    Edalat(tc,pc,ω)

# Edalat relation for vapor pressure

Calculates vapor pressure of a fluid at arbitrary temperatures using a CSP relationship by [1]. Requires a chemical’s critical temperature, pressure, and acentric factor. Claimed to have a higher accuracy than the Lee-Kesler relationship.

The vapor pressure is given by:

    ln(Pᵣ) = f⁰ + f¹ω + f²ω^2
    f⁰ = 6.83377 + -5.76051/Tᵣ + 0.90654*log(Tᵣ) + -1.16906*Tᵣ^1.9
    f¹ = 5.32034 + -28.1460/Tᵣ + -58.0352*log(Tᵣ) + 23.57466*Tᵣ^1.9
    f² = 18.19967 + 16.33839/Tᵣ + 65.6995*log(Tᵣ) + -35.9739*Tᵣ^1.9


## References

(1, 2) Edalat, M., R. B. Bozar-Jomehri, and G. A. Mansoori. “Generalized Equation Predicts Vapor Pressure of Hydrocarbons.” Oil and Gas Journal; 91:5 (February 1, 1993).
"""

struct Edalat{T} <: SaturationModel
    tc::T
    pc::T
    ω::T
end

function Edalat(;tc,pc,ω)
    _tc = ThermoState.normalize_units(tc)
    _pc = ThermoState.normalize_units(pc)
    return Edalat(_tc,_pc,ω)
end

function Edalat(model::ThermoModel)
    _tc = temperature(model,CriticalPoint(),u"K")
    _pc = pressure(model,CriticalPoint(),u"Pa")
    _ω = acentric_factor(model)
    return Edalat(_tc,_pc,_ω)
end

function temperature(model::Edalat,st::CriticalPoint,unit=u"K")
    return convert_unit(u"K",unit,model.tc)
end

function pressure(model::Edalat,st::CriticalPoint,unit=u"Pa")
    return convert_unit(u"Pa",unit,model.pc)
end

function acentric_factor(model::Edalat)
    return model.ω
end

function pressure_impl(mt::SingleSatT,model::Edalat,_t)
    t = normalize_units(_t)
    Pc = model.pc
    τ = 1 - t/model.tc
    ω = model.ω
    τ15 = sqrt(τ)*τ
    τ25 = τ15*τ
    τ5 = τ25*τ25
    tau = 1. - T/Tc
    a = -6.1559 - 4.0855ω
    c = -0.8747 - 7.8874ω
    d = 1/(-0.4893 - 0.9912ω + 3.1551ω^2)
    b = 1.5737 - 1.0540ω - 4.4365E-3*d
    τ15 = sqrt(τ)*τ
    τ3 = τ15*τ15
    lnPr = (a*τ + b*τ15 + c*τ3 + d*τ3*τ3)*model.tc/t
    return exp(lnPr)*Pc
end

initial_temperature(model::Edalat,p) =  critical_sat_interpolation(model,p)


"""
    WaterSat()
    Antoine(a,b,c;from =u"K",to = u"Pa",base=ℯ)

## Water saturation acording to IAPWS-IF97

Calculates vapor pressure of water acording IAPWS Industrial Formulation 

The implementation is from https://github.com/braamvandyk/SteamTables.jl
"""
struct WaterSat <: SaturationModel end

const watersat_data = (;n = [0.116_705_214_527_67E4,
    -0.724_213_167_032_06E6,
    -0.170_738_469_400_92E2,
    0.120_208_247_024_70E5,
    -0.323_255_503_223_33E7,
    0.149_151_086_135_30E2,
    -0.482_326_573_615_91E4,
    0.405_113_405_420_57E6,
    -0.238_555_575_678_49,
    0.650_175_348_447_98E3]
)

function pressure_impl(mt::SingleSatT,::WaterSat,t)
    n = watersat_data.n
    Θ = t + n[9]/(t - n[10])
    A =      Θ^2 + n[1]*Θ + n[2]
    B = n[3]*Θ^2 + n[4]*Θ + n[5]
    C = n[6]*Θ^2 + n[7]*Θ + n[8]
    P = (2C / (-B + √(B^2 - 4*A*C)))^4
    return convert_unit(u"MPa",u"Pa",P)
end

function temperature_impl(mt::SingleSatP,::WaterSat,p)
    n = watersat_data.n
    P = convert_unit(u"Pa",u"MPa",p)
    β = P^0.25
    E =      β^2 + n[3]*β + n[6]
    F = n[1]*β^2 + n[4]*β + n[7]
    G = n[2]*β^2 + n[5]*β + n[8]
    D = 2G / (-F - √(F^2 - 4*E*G))
    T = (n[10]+D-√((n[10]+D)^2 - 4(n[9]+n[10]*D)))/2.0
    return T
end




 
