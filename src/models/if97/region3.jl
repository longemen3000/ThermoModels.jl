
const Region3_vt =(
    n = [0.106_580_700_285_13E1,
    -0.157_328_452_902_39E2,
     0.209_443_969_743_07E2,
    -0.768_677_078_787_16E1,
     0.261_859_477_879_54E1,
    -0.280_807_811_486_20E1,
     0.120_533_696_965_17E1,
    -0.845_668_128_125_02E-2,
    -0.126_543_154_777_14E1,
    -0.115_244_078_066_81E1,
     0.885_210_439_843_18,
    -0.642_077_651_816_07,
     0.384_934_601_866_71,
    -0.852_147_088_242_06,
     0.489_722_815_418_77E1,
    -0.305_026_172_569_65E1,
     0.394_205_368_791_54E-1,
     0.125_584_084_243_08,
    -0.279_993_296_987_10,
     0.138_997_995_694_60E1,
    -0.201_899_150_235_70E1,
    -0.821_476_371_739_63E-2,
    -0.475_960_357_349_23,
     0.439_840_744_735_00E-1,
    -0.444_764_354_287_39,
     0.905_720_707_197_33,
     0.705_224_500_879_67,
     0.107_705_126_263_32,
    -0.329_136_232_589_54,
    -0.508_710_620_411_58,
    -0.221_754_008_730_96E-1,
     0.942_607_516_650_92E-1,
     0.164_362_784_479_61,
    -0.135_033_722_413_48E-1,
    -0.148_343_453_524_72E-1,
     0.579_229_536_280_84E-3,
     0.323_089_047_037_11E-2,
     0.809_648_029_962_15E-4,
    -0.165_576_797_950_37E-3,
    -0.449_238_990_618_15E-4],

I = [0,
     0,
     0,
     0,
     0,
     0,
     0,
     0,
     1,
     1,
     1,
     1,
     2,
     2,
     2,
     2,
     2,
     2,
     3,
     3,
     3,
     3,
     3,
     4,
     4,
     4,
     4,
     5,
     5,
     5,
     6,
     6,
     6,
     7,
     8,
     9,
     9,
     10,
     10,
     11],

J = [0,
     0,
     1,
     2,
     7,
     10,
     12,
     23,
     2,
     6,
     15,
     17,
     0,
     2,
     6,
     7,
     22,
     26,
     0,
     2,
     4,
     16,
     26,
     0,
     2,
     4,
     26,
     1,
     3,
     26,
     0,
     2,
     26,
     2,
     26,
     2,
     26,
     0,
     1,
     26])

     #==
function Region3_ρ(Output::Symbol, _ρ, _T)
    ρ = normalize_units(_ρ)
    T = normalize_units(_T)
    ρstar = IF97_ρc      #kg/m3
    Tstar = IF97_Tc   #K
    n = Region3_vt.n
    I = Region3_vt.I 
    J = Region3_vt.J

    δ = ρ / ρstar
    τ = Tstar / T

    ϕ    =  n[1]*log(δ) + sum(n[i]*(δ^I[i])*(τ^J[i]) for i=2:40)
    ϕ_δ  =  n[1]/δ + sum(n[i]*I[i]*(δ^(I[i]-1))*(τ^J[i]) for i=2:40)
    ϕ_δδ = -n[1]/(δ^2) + sum(n[i]*I[i]*(I[i]-1)*(δ^(I[i]-2))*(τ^J[i]) for i=2:40)
    ϕ_τ  =  sum(n[i]*(δ^I[i])*J[i]*(τ^(J[i]-1)) for i=2:40)
    ϕ_ττ =  sum(n[i]*(δ^I[i])*J[i]*(J[i]-1)*(τ^(J[i]-2)) for i=2:40)
    ϕ_δτ =  sum(n[i]*I[i]*(δ^(I[i]-1))*J[i]*(τ^(J[i]-1)) for i=2:40)

    if Output == :SpecificG            #kJ/kg
        return IF97_R*T*(δ*ϕ_δ + ϕ)
    elseif Output == :SpecificF        #kJ/kg
        return IF97_R*T*ϕ
    elseif Output == :Pressure         #MPa
        return ρ*IF97_R*T*δ*ϕ_δ/1000
    elseif Output == :SpecificU        #kJ/kg
        return IF97_R*T*τ*ϕ_τ
    elseif Output == :SpecificS        #kJ/kgK
        return IF97_R*(τ*ϕ_τ - ϕ)
    elseif Output == :SpecificH        #kJ/kg
        return IF97_R*T*(τ*ϕ_τ + δ*ϕ_δ)
    elseif Output == :SpecificCP       #kJ/kgK
        return IF97_R*(-τ^2*ϕ_ττ+((δ*ϕ_δ - δ*τ*ϕ_δτ)^2)/(2*δ*ϕ_δ + δ^2*ϕ_δδ))
    elseif Output == :SpecificCV       #kJ/kgK
        return IF97_R*(-τ^2*ϕ_ττ)
    elseif Output == :SpeedOfSound     #m/s
        return sqrt(1000*IF97_R*T*(2*δ*ϕ_δ + δ^2*ϕ_δδ - ((δ*ϕ_δ - δ*τ*ϕ_δτ)^2)/(τ^2*ϕ_ττ)))
    else
        throw(DomainError(Output, "Unknown value requested."))
    end
end
==#
function mass_gibbs_impl(mt::SingleVT,model::IF97Region{:r3}, _ρ, _T)
    ρ = normalize_units(_ρ)
    T = normalize_units(_T)
    ρstar = IF97_ρc      #kg/m3
    Tstar = IF97_Tc   #K
    n = Region3_vt.n
    I = Region3_vt.I 
    J = Region3_vt.J
    δ = ρ / ρstar
    τ = Tstar / T
    ϕ    =  n[1]*log(δ) + sum(n[i]*(δ^I[i])*(τ^J[i]) for i=2:40)
    ϕ_δ  =  n[1]/δ + sum(n[i]*I[i]*(δ^(I[i]-1))*(τ^J[i]) for i=2:40)
    res = IF97_R*T*(δ*ϕ_δ + ϕ)
    return convert_unit(u"kJ/kg",u"J/kg",res)
end

function mass_helmholtz_impl(mt::SingleVT,model::IF97Region{:r3}, _ρ, _T)
    ρ = normalize_units(_ρ)
    T = normalize_units(_T)
    ρstar = IF97_ρc      #kg/m3
    Tstar = IF97_Tc   #K
    n = Region3_vt.n
    I = Region3_vt.I 
    J = Region3_vt.J
    δ = ρ / ρstar
    τ = Tstar / T
    ϕ    =  n[1]*log(δ) + sum(n[i]*(δ^I[i])*(τ^J[i]) for i=2:40)

    res = IF97_R*T*ϕ
    return convert_unit(u"kJ/kg",u"J/kg",res)
end

function pressure_impl(mt::SingleVT,model::IF97Region{:r3}, _ρ, _T)
    ρ = normalize_units(_ρ)
    T = normalize_units(_T)
    ρstar = IF97_ρc      #kg/m3
    Tstar = IF97_Tc   #K
    n = Region3_vt.n
    I = Region3_vt.I 
    J = Region3_vt.J
    δ = ρ / ρstar
    τ = Tstar / T
    ϕ_δ  =  n[1]/δ + sum(n[i]*I[i]*(δ^(I[i]-1))*(τ^J[i]) for i=2:40)
    res = ρ*IF97_R*T*δ*ϕ_δ/1000
    return convert_unit(u"MPa",u"Pa",res)
end

function mass_internal_energy_impl(mt::SingleVT,model::IF97Region{:r3}, _ρ, _T)
    ρ = normalize_units(_ρ)
    T = normalize_units(_T)
    ρstar = IF97_ρc      #kg/m3
    Tstar = IF97_Tc   #K
    n = Region3_vt.n
    I = Region3_vt.I 
    J = Region3_vt.J
    δ = ρ / ρstar
    τ = Tstar / T
    ϕ    =  n[1]*log(δ) + sum(n[i]*(δ^I[i])*(τ^J[i]) for i=2:40)
    ϕ_τ  =  sum(n[i]*(δ^I[i])*J[i]*(τ^(J[i]-1)) for i=2:40)
    res = IF97_R*T*τ*ϕ_τ
    return convert_unit(u"kJ/kg",u"J/kg",res)
end

function mass_entropy_impl(mt::SingleVT,model::IF97Region{:r3}, _ρ, _T)
    ρ = normalize_units(_ρ)
    T = normalize_units(_T)
    ρstar = IF97_ρc      #kg/m3
    Tstar = IF97_Tc   #K
    n = Region3_vt.n
    I = Region3_vt.I 
    J = Region3_vt.J
    δ = ρ / ρstar
    τ = Tstar / T
    ϕ    =  n[1]*log(δ) + sum(n[i]*(δ^I[i])*(τ^J[i]) for i=2:40)
    ϕ_τ  =  sum(n[i]*(δ^I[i])*J[i]*(τ^(J[i]-1)) for i=2:40)
    res = IF97_R*(τ*ϕ_τ - ϕ)
    return convert_unit(u"kJ/(kg*K)",u"J/(kg*K)",res)
end

function mass_enthalpy_impl(mt::SingleVT,model::IF97Region{:r3}, _ρ, _T)
    ρ = normalize_units(_ρ)
    T = normalize_units(_T)
    ρstar = IF97_ρc      #kg/m3
    Tstar = IF97_Tc   #K
    n = Region3_vt.n
    I = Region3_vt.I 
    J = Region3_vt.J
    δ = ρ / ρstar
    τ = Tstar / T
    ϕ    =  n[1]*log(δ) + sum(n[i]*(δ^I[i])*(τ^J[i]) for i=2:40)
    ϕ_δ  =  n[1]/δ + sum(n[i]*I[i]*(δ^(I[i]-1))*(τ^J[i]) for i=2:40)
    ϕ_τ  =  sum(n[i]*(δ^I[i])*J[i]*(τ^(J[i]-1)) for i=2:40)
    res = IF97_R*T*(τ*ϕ_τ + δ*ϕ_δ)
    return convert_unit(u"kJ/kg",u"J/kg",res)
end

function mass_cp_impl(mt::SingleVT,model::IF97Region{:r3}, _ρ, _T)
    ρ = normalize_units(_ρ)
    T = normalize_units(_T)
    ρstar = IF97_ρc      #kg/m3
    Tstar = IF97_Tc   #K
    n = Region3_vt.n
    I = Region3_vt.I 
    J = Region3_vt.J
    δ = ρ / ρstar
    τ = Tstar / T
    ϕ_ττ =  sum(n[i]*(δ^I[i])*J[i]*(J[i]-1)*(τ^(J[i]-2)) for i=2:40)
    ϕ_δ  =  n[1]/δ + sum(n[i]*I[i]*(δ^(I[i]-1))*(τ^J[i]) for i=2:40)
    ϕ_δτ =  sum(n[i]*I[i]*(δ^(I[i]-1))*J[i]*(τ^(J[i]-1)) for i=2:40)
    ϕ_δδ = -n[1]/(δ^2) + sum(n[i]*I[i]*(I[i]-1)*(δ^(I[i]-2))*(τ^J[i]) for i=2:40)

    res = IF97_R*(-τ^2*ϕ_ττ+((δ*ϕ_δ - δ*τ*ϕ_δτ)^2)/(2*δ*ϕ_δ + δ^2*ϕ_δδ))
    return convert_unit(u"kJ/(kg*K)",u"J/(kg*K)",res)
end

function mass_cv_impl(mt::SingleVT,model::IF97Region{:r3}, _ρ, _T)
    ρ = normalize_units(_ρ)
    T = normalize_units(_T)
    ρstar = IF97_ρc      #kg/m3
    Tstar = IF97_Tc   #K
    n = Region3_vt.n
    I = Region3_vt.I 
    J = Region3_vt.J
    δ = ρ / ρstar
    τ = Tstar / T
    ϕ_ττ =  sum(n[i]*(δ^I[i])*J[i]*(J[i]-1)*(τ^(J[i]-2)) for i=2:40)

    res = IF97_R*(-τ^2*ϕ_ττ)
    return convert_unit(u"kJ/(kg*K)",u"J/(kg*K)",res)
end

function sound_speed_impl(mt::SingleVT,model::IF97Region{:r3}, _ρ, _T)
    ρ = normalize_units(_ρ)
    T = normalize_units(_T)
    ρstar = IF97_ρc      #kg/m3
    Tstar = IF97_Tc   #K
    n = Region3_vt.n
    I = Region3_vt.I 
    J = Region3_vt.J
    δ = ρ / ρstar
    τ = Tstar / T
    ϕ_ττ =  sum(n[i]*(δ^I[i])*J[i]*(J[i]-1)*(τ^(J[i]-2)) for i=2:40)
    ϕ_δ  =  n[1]/δ + sum(n[i]*I[i]*(δ^(I[i]-1))*(τ^J[i]) for i=2:40)
    ϕ_δτ =  sum(n[i]*I[i]*(δ^(I[i]-1))*J[i]*(τ^(J[i]-1)) for i=2:40)
    ϕ_δδ = -n[1]/(δ^2) + sum(n[i]*I[i]*(I[i]-1)*(δ^(I[i]-2))*(τ^J[i]) for i=2:40)

    res = sqrt(1000*IF97_R*T*(2*δ*ϕ_δ + δ^2*ϕ_δδ - ((δ*ϕ_δ - δ*τ*ϕ_δτ)^2)/(τ^2*ϕ_ττ)))
    return res
end

#==
function Region3(Output::Symbol, P, T)
    ρ0 = 1000*P/(R*T) #Starting value from ideal gas

    f(ρ) = Region3_ρ(:Pressure, ρ, T) - P
    ρ = Roots.find_zero(f, ρ0)

    if Output == :SpecificV
        return 1.0 / ρ
    else
        return Region3_ρ(Output, ρ, T)
    end
end
==#

for op in [:mass_gibbs_impl,:mass_helmholtz_impl,
    :mass_enthalpy_impl,:mass_internal_energy_impl,:mass_entropy_impl,
    :mass_cv_impl,:mass_cp_impl,:sound_speed_impl]

    @eval begin
        function $op(mt::SinglePT,model::IF97Region{:r3},p,t)
            P =convert_unit(u"Pa",u"MPa",p)
            T = normalize_units(t)
            ρ0 = 1000*P/(IF97_R*T) #Starting value from ideal gas
            f(ρ) = pressure_impl(QuickStates.ρt(), model,ρ, T) - normalize_units(p)
            ρ = Roots.find_zero(f, ρ0)
            return $op(QuickStates.ρt(), model,ρ, T)
        end
    end

end

function mass_volume_impl(mt::SinglePT,model::IF97Region{:r3},p,t)
    P =convert_unit(u"Pa",u"MPa",p)
    T = normalize_units(t)
    ρ0 = 1000*P/(IF97_R*T) #Starting value from ideal gas
    f(ρ) = pressure_impl(QuickStates.ρt(), model,ρ, T) - normalize_units(p)
    ρ = Roots.find_zero(f, ρ0)
    return inv(ρ)
end

