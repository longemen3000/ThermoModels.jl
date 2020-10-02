const Region1_pt = (
n = [0.146_329_712_131_67,
-0.845_481_871_691_14,
-0.375_636_036_720_40E1,
 0.338_551_691_683_85E1,
-0.957_919_633_878_72,
 0.157_720_385_132_28,
-0.166_164_171_995_01E-1,
 0.812_146_299_835_68E-3,
 0.283_190_801_238_04E-3,
-0.607_063_015_658_74E-3,
-0.189_900_682_184_19E-1,
-0.325_297_487_705_05E-1,
-0.218_417_171_754_14E-1,
-0.528_383_579_699_30E-4,
-0.471_843_210_732_67E-3,
-0.300_017_807_930_26E-3,
 0.476_613_939_069_87E-4,
-0.441_418_453_308_46E-5,
-0.726_949_962_975_94E-15,
-0.316_796_448_450_54E-4,
-0.282_707_979_853_12E-5,
-0.852_051_281_201_03E-9,
-0.224_252_819_080_00E-5,
-0.651_712_228_956_01E-6,
-0.143_417_299_379_24E-12,
-0.405_169_968_601_17E-6,
-0.127_343_017_416_41E-8,
-0.174_248_712_306_34E-9,
-0.687_621_312_955_31E-18,
 0.144_783_078_285_21E-19,
 0.263_357_816_627_95E-22,
-0.119_476_226_400_71E-22,
 0.182_280_945_814_04E-23,
-0.935_370_872_924_58E-25],

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
 1,
 1,
 2,
 2,
 2,
 2,
 2,
 3,
 3,
 3,
 4,
 4,
 4,
 5,
 8,
 8,
 21,
 23,
 29,
 30,
 31,
 32],

J = [-2,
 -1,
 0,
 1,
 2,
 3,
 4,
 5,
 -9,
 -7,
 -1,
 0,
 1,
 3,
 -3,
 0,
 1,
 3,
 17,
 -4,
 0,
 6,
 -5,
 -2,
 10,
 -8,
 -11,
 -6,
 -29,
 -31,
 -38,
 -39,
 -40,
 -41])
"""
    IF97Region{:r1} <: ThermoModel

Returns all the property values in region 1.
273.15K ≤ T ≤ 623.15K  Psat(T) ≤ P ≤ 100MPa
Pressures in MPa and temperature in [K]
"""
const IF97Region1 =  IF97Region{:r1} 

function mass_gibbs_impl(mt::SinglePT,model::IF97Region{:r1},p,t)
    P = convert_unit(u"Pa",u"MPa",p)
    T = normalize_units(t)
    Pstar = 16.53   #MPa
    Tstar = 1386.0  #K
    n = Region1_pt.n
    I = Region1_pt.I
    J = Region1_pt.J
    π = P / Pstar
    τ = Tstar / T
    ππ = 7.1 - π
    ττ = τ - 1.222
    γ    = sum(n[i]*(ππ^I[i])*(ττ^J[i]) for i=1:34)
    res =  IF97_R*T*γ
    return convert_unit(u"kJ/kg",u"J/kg",res)
end

function mass_helmholtz_impl(mt::SinglePT,model::IF97Region{:r1},p,t)
    P = convert_unit(u"Pa",u"MPa",p)
    T = normalize_units(t)
    Pstar = 16.53   #MPa
    Tstar = 1386.0  #K
    n = Region1_pt.n
    I = Region1_pt.I
    J = Region1_pt.J
    π = P / Pstar
    τ = Tstar / T
    ππ = 7.1 - π
    ττ = τ - 1.222
    γ    = sum(n[i]*(ππ^I[i])*(ττ^J[i]) for i=1:34)
    γ_π  = sum(-n[i]*I[i]*(ππ^(I[i]-1))*(ττ^J[i]) for i=1:34)

    res =  IF97_R*T*(γ - π*γ_π/1000)
    return convert_unit(u"kJ/kg",u"J/kg",res)
end

function mass_volume_impl(mt::SinglePT,model::IF97Region{:r1},p,t)
    P = convert_unit(u"Pa",u"MPa",p)
    T = normalize_units(t)
    Pstar = 16.53   #MPa
    Tstar = 1386.0  #K
    n = Region1_pt.n
    I = Region1_pt.I
    J = Region1_pt.J
    π = P / Pstar
    τ = Tstar / T
    ππ = 7.1 - π
    ττ = τ - 1.222
    γ_π  = sum(-n[i]*I[i]*(ππ^(I[i]-1))*(ττ^J[i]) for i=1:34)

    res = IF97_R*T*π*γ_π/P/1000
    return res
end

function mass_internal_energy_impl(mt::SinglePT,model::IF97Region{:r1},p,t)
    P = convert_unit(u"Pa",u"MPa",p)
    T = normalize_units(t)
    Pstar = 16.53   #MPa
    Tstar = 1386.0  #K
    n = Region1_pt.n
    I = Region1_pt.I
    J = Region1_pt.J
    π = P / Pstar
    τ = Tstar / T
    ππ = 7.1 - π
    ττ = τ - 1.222
    γ_π  = sum(-n[i]*I[i]*(ππ^(I[i]-1))*(ττ^J[i]) for i=1:34)
    γ_τ  = sum(n[i]*(ππ^I[i])*J[i]*(ττ^(J[i]-1)) for i=1:34)
    res =  IF97_R*T*(τ*γ_τ - π*γ_π)
    return convert_unit(u"kJ/kg",u"J/kg",res)
end

function mass_entropy_impl(mt::SinglePT,model::IF97Region{:r1},p,t)
    P = convert_unit(u"Pa",u"MPa",p)
    T = normalize_units(t)
    Pstar = 16.53   #MPa
    Tstar = 1386.0  #K
    n = Region1_pt.n
    I = Region1_pt.I
    J = Region1_pt.J
    π = P / Pstar
    τ = Tstar / T
    ππ = 7.1 - π
    ττ = τ - 1.222
    γ    = sum(n[i]*(ππ^I[i])*(ττ^J[i]) for i=1:34)
    γ_τ  = sum(n[i]*(ππ^I[i])*J[i]*(ττ^(J[i]-1)) for i=1:34)
    res =  IF97_R*(τ*γ_τ - γ)
    return convert_unit(u"kJ/(kg*K)",u"J/(kg*K)",res)
end

function mass_enthalpy_impl(mt::SinglePT,model::IF97Region{:r1},p,t)
    P = convert_unit(u"Pa",u"MPa",p)
    T = normalize_units(t)
    Pstar = 16.53   #MPa
    Tstar = 1386.0  #K
    n = Region1_pt.n
    I = Region1_pt.I
    J = Region1_pt.J
    π = P / Pstar
    τ = Tstar / T
    ππ = 7.1 - π
    ττ = τ - 1.222
    γ    = sum(n[i]*(ππ^I[i])*(ττ^J[i]) for i=1:34)
    γ_τ  = sum(n[i]*(ππ^I[i])*J[i]*(ττ^(J[i]-1)) for i=1:34)

    res = IF97_R*T*τ*γ_τ
    return convert_unit(u"kJ/kg",u"J/kg",res)
end

function mass_cp_impl(mt::SinglePT,model::IF97Region{:r1},p,t)
    P = convert_unit(u"Pa",u"MPa",p)
    T = normalize_units(t)
    Pstar = 16.53   #MPa
    Tstar = 1386.0  #K
    n = Region1_pt.n
    I = Region1_pt.I
    J = Region1_pt.J
    π = P / Pstar
    τ = Tstar / T
    ππ = 7.1 - π
    ττ = τ - 1.222
    γ_ττ = sum(n[i]*(ππ^I[i])*J[i]*(J[i]-1)*(ττ^(J[i]-2)) for i=1:34)

    res = return -IF97_R*τ^2*γ_ττ
    return convert_unit(u"kJ/(kg*K)",u"J/(kg*K)",res)

end

function mass_cv_impl(mt::SinglePT,model::IF97Region{:r1},p,t)
    P = convert_unit(u"Pa",u"MPa",p)
    T = normalize_units(t)
    Pstar = 16.53   #MPa
    Tstar = 1386.0  #K
    n = Region1_pt.n
    I = Region1_pt.I
    J = Region1_pt.J
    π = P / Pstar
    τ = Tstar / T
    ππ = 7.1 - π
    ττ = τ - 1.222
    γ_ππ = sum(n[i]*I[i]*(I[i]-1)*(ππ^(I[i]-2))*(ττ^J[i]) for i=1:34)
    γ_πτ = sum(-n[i]*I[i]*(ππ^(I[i]-1))*J[i]*(ττ^(J[i]-1)) for i=1:34)
    γ_π  = sum(-n[i]*I[i]*(ππ^(I[i]-1))*(ττ^J[i]) for i=1:34)
    γ_ττ = sum(n[i]*(ππ^I[i])*J[i]*(J[i]-1)*(ττ^(J[i]-2)) for i=1:34)

    res =  IF97_R*((-τ^2)*γ_ττ + (γ_π - τ*γ_πτ)^2/γ_ππ)
    return convert_unit(u"kJ/(kg*K)",u"J/(kg*K)",res)

end

function sound_speed_impl(mt::SinglePT,model::IF97Region{:r1},p,t)
    P = convert_unit(u"Pa",u"MPa",p)
    T = normalize_units(t)
    Pstar = 16.53   #MPa
    Tstar = 1386.0  #K
    n = Region1_pt.n
    I = Region1_pt.I
    J = Region1_pt.J
    π = P / Pstar
    τ = Tstar / T
    ππ = 7.1 - π
    ττ = τ - 1.222
    γ_ππ = sum(n[i]*I[i]*(I[i]-1)*(ππ^(I[i]-2))*(ττ^J[i]) for i=1:34)
    γ_ττ = sum(n[i]*(ππ^I[i])*J[i]*(J[i]-1)*(ττ^(J[i]-2)) for i=1:34)
    γ_πτ = sum(-n[i]*I[i]*(ππ^(I[i]-1))*J[i]*(ττ^(J[i]-1)) for i=1:34)
    γ_π  = sum(-n[i]*I[i]*(ππ^(I[i]-1))*(ττ^J[i]) for i=1:34)

    res =  return √(1000*IF97_R*T*(γ_π)^2/((γ_π - τ*γ_πτ)^2/τ^2/γ_ττ - γ_ππ))
    return res
end


const Region1_t_ph =(
n = [-0.238_724_899_245_21E3,
0.404_211_886_379_45E3,
0.113_497_468_817_18E3,
-0.584_576_160_480_39E1,
-0.152_854_824_131_40E-3,
-0.108_667_076_953_77E-5,
-0.133_917_448_726_02E2,
0.432_110_391_835_59E2,
-0.540_100_671_705_06E2,
0.305_358_922_039_16E2,
-0.659_647_494_236_38E1,
0.939_654_008_783_63E-2,
0.115_736_475_053_40E-6,
-0.258_586_412_820_73E-4,
-0.406_443_630_847_99E-8,
0.664_561_861_916_35E-7,
0.806_707_341_030_27E-10,
-0.934_777_712_139_47E-12,
0.582_654_420_206_01E-14,
-0.150_201_859_535_03E-16],

I = [0,
0,
0,
0,
0,
0,
1,
1,
1,
1,
1,
1,
1,
2,
2,
3,
3,
4,
5,
6],

J = [0,
1,
2,
6,
22,
32,
0,
1,
2,
3,
4,
10,
32,
10,
32,
10,
32,
32,
32,
32])



function temperature_impl(mt::SinglePH,model::IF97Region{:r1},p,_h)
    P = convert_unit(u"Pa",u"MPa",p)
    h = convert_unit(u"J/kg",u"kJ/kg",normalize_units(_h))
    hstar = 2500 #kJ/kg

    n = Region1_t_ph.n
    I = Region1_t_ph.I
    J = Region1_t_ph.J
    η = h / hstar
    return sum(n[i]*P^I[i]*(η+1)^J[i] for i=1:20)
end


const Region1_t_ps =(
    n = [0.174_782_680_583_07E3,
    0.348_069_308_928_73E2,
    0.652_925_849_784_55E1,
    0.330_399_817_754_89,
   -0.192_813_829_231_96E-6,
   -0.249_091_972_445_73E-22,
   -0.261_076_364_893_32,
    0.225_929_659_815_86,
   -0.642_564_633_952_26E-1,
    0.788_762_892_705_26E-2,
    0.356_721_106_073_66E-9,
    0.173_324_969_948_95E-23,
    0.566_089_006_548_37E-3,
   -0.326_354_831_397_17E-3,
    0.447_782_866_906_32E-4,
   -0.513_221_569_085_07E-9,
   -0.425_226_570_422_07E-25,
    0.264_004_413_606_89E-12,
    0.781_246_004_597_23E-28,
   -0.307_321_999_036_68E-30],

I = [0,
    0,
    0,
    0,
    0,
    0,
    1,
    1,
    1,
    1,
    1,
    1,
    2,
    2,
    2,
    2,
    2,
    3,
    3,
    4],

J = [0,
    1,
    2,
    3,
    11,
    31,
    0,
    1,
    2,
    3,
    12,
    31,
    0,
    1,
    2,
    9,
    31,
    10,
    32,
    32]

)

function temperature_impl(mt::SinglePS,model::IF97Region{:r1},p,_s)
    P = convert_unit(u"Pa",u"MPa",p)
    s = convert_unit(u"J/(kg*K)",u"kJ/(kg*K)",normalize_units(_s))
    n = Region1_t_ps.n
    I = Region1_t_ps.I
    J = Region1_t_ps.J

    return sum(n[i]*P^I[i]*(s+2)^J[i] for i=1:20)
end

