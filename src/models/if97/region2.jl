const Region2_pt = ( no = [-0.969_276_865_002_17E1,
    0.100_866_559_680_18E2,
    -0.560_879_112_830_20E-2,
    0.714_527_380_814_55E-1,
    -0.407_104_982_239_28,
    0.142_408_191_714_44E1,
    -0.438_395_113_194_50E1,
    -0.284_086_324_607_72,
    0.212_684_637_533_07E-1],

    Jo = [0,
    1,
    -5,
    -4,
    -3,
    -2,
    -1,
    2,
    3],

    nr = [-0.177_317_424_732_13E-2,
    -0.178_348_622_923_58E-1,
    -0.459_960_136_963_65E-1,
    -0.575_812_590_834_32E-1,
    -0.503_252_787_279_30E-1,
    -0.330_326_416_702_03E-4,
    -0.189_489_875_163_15E-3,
    -0.393_927_772_433_55E-2,
    -0.437_972_956_505_73E-1,
    -0.266_745_479_140_87E-4,
    0.204_817_376_923_09E-7,
    0.438_706_672_844_35E-6,
    -0.322_776_772_385_70E-4,
    -0.150_339_245_421_48E-2,
    -0.406_682_535_626_49E-1,
    -0.788_473_095_593_67E-9,
    0.127_907_178_522_85E-7,
    0.482_253_727_185_07E-6,
    0.229_220_763_376_61E-5,
    -0.167_147_664_510_61E-10,
    -0.211_714_723_213_55E-2,
    -0.238_957_419_341_04E2,
    -0.590_595_643_242_70E-17,
    -0.126_218_088_991_01E-5,
    -0.389_468_424_357_39E-1,
    0.112_562_113_604_59E-10,
    -0.823_113_408_979_98E1,
    0.198_097_128_020_88E-7,
    0.104_069_652_101_74E-18,
    -0.102_347_470_959_29E-12,
    -0.100_181_793_795_11E-8,
    -0.808_829_086_469_85E-10,
    0.106_930_318_794_09,
    -0.336_622_505_741_71,
    0.891_858_453_554_21E-24,
    0.306_293_168_762_32E-12,
    -0.420_024_676_982_08E-5,
    -0.590_560_296_856_39E-25,
    0.378_269_476_134_57E-5,
    -0.127_686_089_346_81E-14,
    0.730_876_105_950_61E-28,
    0.554_147_153_507_78E-16,
    -0.943_697_072_412_10E-6],

    Ir = [1,
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
    3,
    3,
    4,
    4,
    4,
    5,
    6,
    6,
    6,
    7,
    7,
    7,
    8,
    8,
    9,
    10,
    10,
    10,
    16,
    16,
    18,
    20,
    20,
    20,
    21,
    22,
    23,
    24,
    24,
    24],

    Jr = [0,
    1,
    2,
    3,
    6,
    1,
    2,
    4,
    7,
    36,
    0,
    1,
    3,
    6,
    35,
    1,
    2,
    3,
    7,
    3,
    16,
    35,
    0,
    11,
    25,
    8,
    36,
    13,
    4,
    10,
    14,
    29,
    50,
    57,
    20,
    35,
    48,
    21,
    53,
    39,
    26,
    40,
    58])

"""
    IF97Region2

    Returns all the property values in region 2.
    273.15K ≤ T ≤ 623.15K  0 ≤ P ≤ Psat(T)
    623.15K <  T ≤ 863.15K  0 <  P ≤ P(T) from B23-model
    863.15K <  T ≤ 1073.15K 0 <  P ≤ 100MPa
    Accuracy in the metastable region is reasonable above 10MPa, only
"""
const IF97Region2 = IF97Region{:r2}

function Region2(Output::Symbol, P, T)
    P = convert_unit(u"Pa",u"MPa",p)
    T = normalize_units(t)
    Pstar = 1.0    #MPa
    Tstar = 540.0  #K
    no = Region2_pt.no
    Jo = Region2_pt.Jo
    nr = Region2_pt.nr
    Ir = Region2_pt.Ir
    Jr = Region2_pt.Jr

    π = P / Pstar
    τ = Tstar / T

    γo    = log(π) + sum(no[i]*τ^Jo[i] for i=1:9)
    γo_π  = 1/π
    γo_ππ = -1/(π^2)
    γo_τ  = sum(no[i]*Jo[i]*τ^(Jo[i]-1) for i=1:9)
    γo_ττ = sum(no[i]*Jo[i]*(Jo[i]-1)*τ^(Jo[i]-2) for i=1:9)
    γo_πτ = 0

    γr    = sum(nr[i]*π^Ir[i]*(τ-0.5)^Jr[i] for i = 1:43)
    γr_π  = sum(nr[i]*Ir[i]*π^(Ir[i]-1)*(τ-0.5)^Jr[i] for i=1:43)
    γr_ππ = sum(nr[i]*Ir[i]*(Ir[i]-1)*π^(Ir[i]-2)*(τ-0.5)^Jr[i] for i=1:43)
    γr_τ  = sum(nr[i]*π^Ir[i]*Jr[i]*(τ-0.5)^(Jr[i]-1) for i=1:43)
    γr_ττ = sum(nr[i]*π^Ir[i]*Jr[i]*(Jr[i]-1)*(τ-0.5)^(Jr[i]-2) for i=1:43)
    γr_πτ = sum(nr[i]*Ir[i]*π^(Ir[i]-1)*Jr[i]*(τ-0.5)^(Jr[i]-1) for i = 1:43)


    if Output == :SpecificG            #kJ/kg
        return IF97_R*T*(γo + γr)
    elseif Output == :SpecificF        #kJ/kg
        return IF97_R*T*(γo + γr - π*(γo_π + γr_π)/1000)
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
        return IF97_R*(-τ^2*(γo_ττ+γr_ττ) - (1+π*γr_π-τ*π*γr_πτ)^2/(1-π^2*γr_ππ))
    elseif Output == :SpeedOfSound     #m/s
        return sqrt((1000*IF97_R*T)*(1+2*π*γr_π+π^2*γr_π^2)/((1-π^2*γr_ππ)+(1+π*γr_π-τ*π*γr_πτ)^2/τ^2/(γo_ττ + γr_ττ)))
    else
        throw(DomainError(Output, "Unknown value requested."))
    end
end


function mass_gibbs_impl(mt::SinglePT,model::IF97Region{:r2},p,t)
    P = convert_unit(u"Pa",u"MPa",p)
    T = normalize_units(t)
    Pstar = 1.0    #MPa
    Tstar = 540.0  #K
    no = Region2_pt.no
    Jo = Region2_pt.Jo
    nr = Region2_pt.nr
    Ir = Region2_pt.Ir
    Jr = Region2_pt.Jr

    π = P / Pstar
    τ = Tstar / T
    γo = log(π) + sum(no[i]*τ^Jo[i] for i=1:9)
    γr = sum(nr[i]*π^Ir[i]*(τ-0.5)^Jr[i] for i = 1:43)
    res = IF97_R*T*(γo + γr)
    return convert_unit(u"kJ/kg",u"J/kg",res)
end

function mass_helmholtz_impl(mt::SinglePT,model::IF97Region{:r2},p,t)
    P = convert_unit(u"Pa",u"MPa",p)
    T = normalize_units(t)
    Pstar = 1.0    #MPa
    Tstar = 540.0  #K
    no = Region2_pt.no
    Jo = Region2_pt.Jo
    nr = Region2_pt.nr
    Ir = Region2_pt.Ir
    Jr = Region2_pt.Jr
    π = P / Pstar
    τ = Tstar / T
    γo    = log(π) + sum(no[i]*τ^Jo[i] for i=1:9)
    γr    = sum(nr[i]*π^Ir[i]*(τ-0.5)^Jr[i] for i = 1:43)
    γo_π  = 1/π
    γr_π  = sum(nr[i]*Ir[i]*π^(Ir[i]-1)*(τ-0.5)^Jr[i] for i=1:43)
    res =  IF97_R*T*(γo + γr - π*(γo_π + γr_π)/1000)
    return convert_unit(u"kJ/kg",u"J/kg",res)
end

function mass_volume_impl(mt::SinglePT,model::IF97Region{:r2},p,t)
    P = convert_unit(u"Pa",u"MPa",p)
    T = normalize_units(t)
    Pstar = 1.0    #MPa
    Tstar = 540.0  #K
    no = Region2_pt.no
    Jo = Region2_pt.Jo
    nr = Region2_pt.nr
    Ir = Region2_pt.Ir
    Jr = Region2_pt.Jr

    π = P / Pstar
    τ = Tstar / T

    γo_π  = 1/π
    γr_π  = sum(nr[i]*Ir[i]*π^(Ir[i]-1)*(τ-0.5)^Jr[i] for i=1:43)

    res = IF97_R*T*π*(γo_π + γr_π)/P/1000
    return res
end

function mass_internal_energy_impl(mt::SinglePT,model::IF97Region{:r2},p,t)
    P = convert_unit(u"Pa",u"MPa",p)
    T = normalize_units(t)
    Pstar = 1.0    #MPa
    Tstar = 540.0  #K
    no = Region2_pt.no
    Jo = Region2_pt.Jo
    nr = Region2_pt.nr
    Ir = Region2_pt.Ir
    Jr = Region2_pt.Jr

    π = P / Pstar
    τ = Tstar / T
    γo_τ  = sum(no[i]*Jo[i]*τ^(Jo[i]-1) for i=1:9)
    γr_τ  = sum(nr[i]*π^Ir[i]*Jr[i]*(τ-0.5)^(Jr[i]-1) for i=1:43)
    γo_π  = 1/π
    γr_π  = sum(nr[i]*Ir[i]*π^(Ir[i]-1)*(τ-0.5)^Jr[i] for i=1:43)

    res = IF97_R*T*(τ*(γo_τ + γr_τ) - π*(γo_π + γr_π))

    return convert_unit(u"kJ/kg",u"J/kg",res)
end

function mass_entropy_impl(mt::SinglePT,model::IF97Region{:r2},p,t)
    P = convert_unit(u"Pa",u"MPa",p)
    T = normalize_units(t)
    Pstar = 1.0    #MPa
    Tstar = 540.0  #K
    no = Region2_pt.no
    Jo = Region2_pt.Jo
    nr = Region2_pt.nr
    Ir = Region2_pt.Ir
    Jr = Region2_pt.Jr

    π = P / Pstar
    τ = Tstar / T
    γo_τ  = sum(no[i]*Jo[i]*τ^(Jo[i]-1) for i=1:9)
    γr_τ  = sum(nr[i]*π^Ir[i]*Jr[i]*(τ-0.5)^(Jr[i]-1) for i=1:43)

    γo    = log(π) + sum(no[i]*τ^Jo[i] for i=1:9)
    γr    = sum(nr[i]*π^Ir[i]*(τ-0.5)^Jr[i] for i = 1:43)

    res = IF97_R*(τ*(γo_τ + γr_τ) - (γo + γr))

    return convert_unit(u"kJ/(kg*K)",u"J/(kg*K)",res)
end

function mass_enthalpy_impl(mt::SinglePT,model::IF97Region{:r2},p,t)
    P = convert_unit(u"Pa",u"MPa",p)
    T = normalize_units(t)
    Pstar = 1.0    #MPa
    Tstar = 540.0  #K
    no = Region2_pt.no
    Jo = Region2_pt.Jo
    nr = Region2_pt.nr
    Ir = Region2_pt.Ir
    Jr = Region2_pt.Jr

    π = P / Pstar
    τ = Tstar / T
    γo_τ  = sum(no[i]*Jo[i]*τ^(Jo[i]-1) for i=1:9)
    γr_τ  = sum(nr[i]*π^Ir[i]*Jr[i]*(τ-0.5)^(Jr[i]-1) for i=1:43)

    res = IF97_R*T*τ*(γo_τ + γr_τ)
    return convert_unit(u"kJ/kg",u"J/kg",res)
end

function mass_cp_impl(mt::SinglePT,model::IF97Region{:r2},p,t)
    P = convert_unit(u"Pa",u"MPa",p)
    T = normalize_units(t)
       Pstar = 1.0    #MPa
    Tstar = 540.0  #K
    no = Region2_pt.no
    Jo = Region2_pt.Jo
    nr = Region2_pt.nr
    Ir = Region2_pt.Ir
    Jr = Region2_pt.Jr

    π = P / Pstar
    τ = Tstar / T
    γo_ττ = sum(no[i]*Jo[i]*(Jo[i]-1)*τ^(Jo[i]-2) for i=1:9)
    γr_ττ = sum(nr[i]*π^Ir[i]*Jr[i]*(Jr[i]-1)*(τ-0.5)^(Jr[i]-2) for i=1:43)

    res = -IF97_R*τ^2*(γo_ττ + γr_ττ)
   
    return convert_unit(u"kJ/(kg*K)",u"J/(kg*K)",res)

end

function mass_cv_impl(mt::SinglePT,model::IF97Region{:r2},p,t)
    P = convert_unit(u"Pa",u"MPa",p)
    T = normalize_units(t)
    Pstar = 1.0    #MPa
    Tstar = 540.0  #K
    no = Region2_pt.no
    Jo = Region2_pt.Jo
    nr = Region2_pt.nr
    Ir = Region2_pt.Ir
    Jr = Region2_pt.Jr

    π = P / Pstar
    τ = Tstar / T
    
    γo_ττ = sum(no[i]*Jo[i]*(Jo[i]-1)*τ^(Jo[i]-2) for i=1:9)
    γr_ττ = sum(nr[i]*π^Ir[i]*Jr[i]*(Jr[i]-1)*(τ-0.5)^(Jr[i]-2) for i=1:43)
    γr_π  = sum(nr[i]*Ir[i]*π^(Ir[i]-1)*(τ-0.5)^Jr[i] for i=1:43)
    γr_πτ = sum(nr[i]*Ir[i]*π^(Ir[i]-1)*Jr[i]*(τ-0.5)^(Jr[i]-1) for i = 1:43)
    γr_ππ = sum(nr[i]*Ir[i]*(Ir[i]-1)*π^(Ir[i]-2)*(τ-0.5)^Jr[i] for i=1:43) 

    res = IF97_R*(-τ^2*(γo_ττ+γr_ττ) - (1+π*γr_π-τ*π*γr_πτ)^2/(1-π^2*γr_ππ))

    return convert_unit(u"kJ/(kg*K)",u"J/(kg*K)",res)
end

function sound_speed_impl(mt::SinglePT,model::IF97Region{:r2},p,t)
    P = convert_unit(u"Pa",u"MPa",p)
    T = normalize_units(t)
    Pstar = 1.0    #MPa
    Tstar = 540.0  #K
    no = Region2_pt.no
    Jo = Region2_pt.Jo
    nr = Region2_pt.nr
    Ir = Region2_pt.Ir
    Jr = Region2_pt.Jr

    π = P / Pstar
    τ = Tstar / T


    γr_π  = sum(nr[i]*Ir[i]*π^(Ir[i]-1)*(τ-0.5)^Jr[i] for i=1:43)
    γr_ππ = sum(nr[i]*Ir[i]*(Ir[i]-1)*π^(Ir[i]-2)*(τ-0.5)^Jr[i] for i=1:43)
    γr_πτ = sum(nr[i]*Ir[i]*π^(Ir[i]-1)*Jr[i]*(τ-0.5)^(Jr[i]-1) for i = 1:43)
    γo_ττ = sum(no[i]*Jo[i]*(Jo[i]-1)*τ^(Jo[i]-2) for i=1:9)
    γr_ττ = sum(nr[i]*π^Ir[i]*Jr[i]*(Jr[i]-1)*(τ-0.5)^(Jr[i]-2) for i=1:43)

    res =  sqrt((1000*IF97_R*T)*(1+2*π*γr_π+π^2*γr_π^2)/
    ((1-π^2*γr_ππ)+(1+π*γr_π-τ*π*γr_πτ)^2/τ^2/(γo_ττ + γr_ττ)))
    return res
end

const Region2m_pt = (
    no = [-0.969_372_683_930_49E1, #Updated from Region_2
           0.100_872_759_700_06E2, #Updated from Region_2
          -0.560_879_112_830_20E-2,
           0.714_527_380_814_55E-1,
          -0.407_104_982_239_28,
           0.142_408_191_714_44E1,
          -0.438_395_113_194_50E1,
          -0.284_086_324_607_72,
           0.212_684_637_533_07E-1],

    Jo = [ 0,
           1,
          -5,
          -4,
          -3,
          -2,
          -1,
           2,
           3],

    nr = [-0.733_622_601_865_06E-2,
          -0.882_238_319_431_46E-1,
          -0.723_345_552_132_45E-1,
          -0.408_131_785_344_55E-2,
           0.200_978_033_802_07E-2,
          -0.530_459_218_986_42E-1,
          -0.761_904_090_869_70E-2,
          -0.634_980_376_573_13E-2,
          -0.860_430_930_285_88E-1,
           0.753_215_815_227_70E-2,
          -0.792_383_754_461_39E-2,
          -0.228_881_607_784_47E-3,
          -0.264_565_014_828_10E-2],

    Ir = [1,
          1,
          1,
          1,
          2,
          2,
          2,
          3,
          3,
          4,
          4,
          5,
          5],

    Jr = [0,
          2,
          5,
          11,
          1,
          7,
          16,
          4,
          16,
          7,
          10,
          9,
          10]

)


function mass_gibbs_impl(mt::SinglePT,model::IF97Region{:r2m},p,t)
    P = convert_unit(u"Pa",u"MPa",p)
    T = normalize_units(t)
    Pstar = 1.0    #MPa
    Tstar = 540.0  #K
    no = Region2m_pt.no
    Jo = Region2m_pt.Jo
    nr = Region2m_pt.nr
    Ir = Region2m_pt.Ir
    Jr = Region2m_pt.Jr

    π = P / Pstar
    τ = Tstar / T
    γo = log(π) + sum(no[i]*τ^Jo[i] for i=1:9)
    γr = sum(nr[i]*π^Ir[i]*(τ-0.5)^Jr[i] for i = 1:13)
    res = IF97_R*T*(γo + γr)
    return convert_unit(u"kJ/kg",u"J/kg",res)
end

function mass_helmholtz_impl(mt::SinglePT,model::IF97Region{:r2m},p,t)
    P = convert_unit(u"Pa",u"MPa",p)
    T = normalize_units(t)
    Pstar = 1.0    #MPa
    Tstar = 540.0  #K
    no = Region2_pt.no
    Jo = Region2_pt.Jo
    nr = Region2_pt.nr
    Ir = Region2_pt.Ir
    Jr = Region2_pt.Jr
    π = P / Pstar
    τ = Tstar / T
    γo    = log(π) + sum(no[i]*τ^Jo[i] for i=1:9)
    γr    = sum(nr[i]*π^Ir[i]*(τ-0.5)^Jr[i] for i = 1:13)
    γo_π  = 1/π
    γr_π  = sum(nr[i]*Ir[i]*π^(Ir[i]-1)*(τ-0.5)^Jr[i] for i=1:13)
    res =  IF97_R*T*(γo + γr - π*(γo_π + γr_π)/1000)
    return convert_unit(u"kJ/kg",u"J/kg",res)
end

function mass_volume_impl(mt::SinglePT,model::IF97Region{:r2m},p,t)
    P = convert_unit(u"Pa",u"MPa",p)
    T = normalize_units(t)
    Pstar = 1.0    #MPa
    Tstar = 540.0  #K
    no = Region2m_pt.no
    Jo = Region2m_pt.Jo
    nr = Region2m_pt.nr
    Ir = Region2m_pt.Ir
    Jr = Region2m_pt.Jr

    π = P / Pstar
    τ = Tstar / T

    γo_π  = 1/π
    γr_π  = sum(nr[i]*Ir[i]*π^(Ir[i]-1)*(τ-0.5)^Jr[i] for i=1:13)

    res = IF97_R*T*π*(γo_π + γr_π)/P/1000
    return res
end

function mass_internal_energy_impl(mt::SinglePT,model::IF97Region{:r2m},p,t)
    P = convert_unit(u"Pa",u"MPa",p)
    T = normalize_units(t)
    Pstar = 1.0    #MPa
    Tstar = 540.0  #K
    no = Region2m_pt.no
    Jo = Region2m_pt.Jo
    nr = Region2m_pt.nr
    Ir = Region2m_pt.Ir
    Jr = Region2m_pt.Jr


    π = P / Pstar
    τ = Tstar / T
    γo_τ  = sum(no[i]*Jo[i]*τ^(Jo[i]-1) for i=1:9)
    γr_τ  = sum(nr[i]*π^Ir[i]*Jr[i]*(τ-0.5)^(Jr[i]-1) for i=1:13)
    γo_π  = 1/π
    γr_π  = sum(nr[i]*Ir[i]*π^(Ir[i]-1)*(τ-0.5)^Jr[i] for i=1:13)

    res = IF97_R*T*(τ*(γo_τ + γr_τ) - π*(γo_π + γr_π))

    return convert_unit(u"kJ/kg",u"J/kg",res)
end

function mass_entropy_impl(mt::SinglePT,model::IF97Region{:r2m},p,t)
    P = convert_unit(u"Pa",u"MPa",p)
    T = normalize_units(t)
    Pstar = 1.0    #MPa
    Tstar = 540.0  #K
    no = Region2m_pt.no
    Jo = Region2m_pt.Jo
    nr = Region2m_pt.nr
    Ir = Region2m_pt.Ir
    Jr = Region2m_pt.Jr


    π = P / Pstar
    τ = Tstar / T
    γo_τ  = sum(no[i]*Jo[i]*τ^(Jo[i]-1) for i=1:9)
    γr_τ  = sum(nr[i]*π^Ir[i]*Jr[i]*(τ-0.5)^(Jr[i]-1) for i=1:13)

    γo    = log(π) + sum(no[i]*τ^Jo[i] for i=1:9)
    γr    = sum(nr[i]*π^Ir[i]*(τ-0.5)^Jr[i] for i = 1:13)

    res = IF97_R*(τ*(γo_τ + γr_τ) - (γo + γr))

    return convert_unit(u"kJ/(kg*K)",u"J/(kg*K)",res)
end

function mass_enthalpy_impl(mt::SinglePT,model::IF97Region{:r2m},p,t)
    P = convert_unit(u"Pa",u"MPa",p)
    T = normalize_units(t)
    Pstar = 1.0    #MPa
    Tstar = 540.0  #K
    no = Region2m_pt.no
    Jo = Region2m_pt.Jo
    nr = Region2m_pt.nr
    Ir = Region2m_pt.Ir
    Jr = Region2m_pt.Jr


    π = P / Pstar
    τ = Tstar / T
    γo_τ  = sum(no[i]*Jo[i]*τ^(Jo[i]-1) for i=1:9)
    γr_τ  = sum(nr[i]*π^Ir[i]*Jr[i]*(τ-0.5)^(Jr[i]-1) for i=1:13)

    res = IF97_R*T*τ*(γo_τ + γr_τ)
    return convert_unit(u"kJ/kg",u"J/kg",res)
end

function mass_cp_impl(mt::SinglePT,model::IF97Region{:r2m},p,t)
    P = convert_unit(u"Pa",u"MPa",p)
    T = normalize_units(t)
       Pstar = 1.0    #MPa
    Tstar = 540.0  #K
    no = Region2m_pt.no
    Jo = Region2m_pt.Jo
    nr = Region2m_pt.nr
    Ir = Region2m_pt.Ir
    Jr = Region2m_pt.Jr


    π = P / Pstar
    τ = Tstar / T
    γo_ττ = sum(no[i]*Jo[i]*(Jo[i]-1)*τ^(Jo[i]-2) for i=1:9)
    γr_ττ = sum(nr[i]*π^Ir[i]*Jr[i]*(Jr[i]-1)*(τ-0.5)^(Jr[i]-2) for i=1:13)

    res = -IF97_R*τ^2*(γo_ττ + γr_ττ)
   
    return convert_unit(u"kJ/(kg*K)",u"J/(kg*K)",res)

end

function mass_cv_impl(mt::SinglePT,model::IF97Region{:r2m},p,t)
    P = convert_unit(u"Pa",u"MPa",p)
    T = normalize_units(t)
    Pstar = 1.0    #MPa
    Tstar = 540.0  #K
    no = Region2m_pt.no
    Jo = Region2m_pt.Jo
    nr = Region2m_pt.nr
    Ir = Region2m_pt.Ir
    Jr = Region2m_pt.Jr


    π = P / Pstar
    τ = Tstar / T
    
    γo_ττ = sum(no[i]*Jo[i]*(Jo[i]-1)*τ^(Jo[i]-2) for i=1:9)
    γr_ττ = sum(nr[i]*π^Ir[i]*Jr[i]*(Jr[i]-1)*(τ-0.5)^(Jr[i]-2) for i=1:13)
    γr_π  = sum(nr[i]*Ir[i]*π^(Ir[i]-1)*(τ-0.5)^Jr[i] for i=1:13)
    γr_πτ = sum(nr[i]*Ir[i]*π^(Ir[i]-1)*Jr[i]*(τ-0.5)^(Jr[i]-1) for i = 1:13)
    γr_ππ = sum(nr[i]*Ir[i]*(Ir[i]-1)*π^(Ir[i]-2)*(τ-0.5)^Jr[i] for i=1:13) 

    res = IF97_R*(-τ^2*(γo_ττ+γr_ττ) - (1+π*γr_π-τ*π*γr_πτ)^2/(1-π^2*γr_ππ))

    return convert_unit(u"kJ/(kg*K)",u"J/(kg*K)",res)
end

function sound_speed_impl(mt::SinglePT,model::IF97Region{:r2m},p,t)
    P = convert_unit(u"Pa",u"MPa",p)
    T = normalize_units(t)
    Pstar = 1.0    #MPa
    Tstar = 540.0  #K
    no = Region2m_pt.no
    Jo = Region2m_pt.Jo
    nr = Region2m_pt.nr
    Ir = Region2m_pt.Ir
    Jr = Region2m_pt.Jr


    π = P / Pstar
    τ = Tstar / T


    γr_π  = sum(nr[i]*Ir[i]*π^(Ir[i]-1)*(τ-0.5)^Jr[i] for i=1:13)
    γr_ππ = sum(nr[i]*Ir[i]*(Ir[i]-1)*π^(Ir[i]-2)*(τ-0.5)^Jr[i] for i=1:13)
    γr_πτ = sum(nr[i]*Ir[i]*π^(Ir[i]-1)*Jr[i]*(τ-0.5)^(Jr[i]-1) for i = 1:13)
    γo_ττ = sum(no[i]*Jo[i]*(Jo[i]-1)*τ^(Jo[i]-2) for i=1:9)
    γr_ττ = sum(nr[i]*π^Ir[i]*Jr[i]*(Jr[i]-1)*(τ-0.5)^(Jr[i]-2) for i=1:13)

    res =  sqrt((1000*IF97_R*T)*(1+2*π*γr_π+π^2*γr_π^2)/
    ((1-π^2*γr_ππ)+(1+π*γr_π-τ*π*γr_πτ)^2/τ^2/(γo_ττ + γr_ττ)))
    return res
end

const Region2a_t_ph = (
    
    n = [0.108_989_523_182_88E4,
         0.849_516_544_955_35E3,
        -0.107_817_480_918_26E3,
         0.331_536_548_012_63E2,
        -0.742_320_167_902_48E1,
         0.117_650_487_243_56E2,
         0.184_457_493_557_90E1,
        -0.417_927_005_496_24E1,
         0.624_781_969_358_12E1,
        -0.173_445_631_081_14E2,
        -0.200_581_768_620_96E3,
         0.271_960_654_737_96E3,
        -0.455_113_182_858_18E3,
         0.309_196_886_047_55E4,
         0.252_266_403_578_72E6,
        -0.617_074_228_683_39E-2,
        -0.310_780_466_295_83,
         0.116_708_730_771_07E2,
         0.128_127_984_040_46E9,
        -0.985_549_096_232_76E9,
         0.282_245_469_730_02E10,
        -0.359_489_714_107_03E10,
         0.172_273_499_131_97E10,
        -0.135_513_342_407_75E5,
         0.128_487_346_646_50E8,
         0.138_657_242_832_26E1,
         0.235_988_325_565_14E6,
        -0.131_052_365_450_54E8,
         0.739_998_354_747_66E4,
        -0.551_966_970_300_60E6,
         0.371_540_859_962_33E7,
         0.191_277_292_396_60E5,
        -0.415_351_648_356_34E6,
        -0.624_598_551_925_07E2],

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
         1,
         1,
         2,
         2,
         2,
         2,
         2,
         2,
         2,
         2,
         3,
         3,
         4,
         4,
         4,
         5,
         5,
         5,
         6,
         6,
         7],

    J = [0,
         1,
         2,
         3,
         7,
         20,
         0,
         1,
         2,
         3,
         7,
         9,
         11,
         18,
         44,
         0,
         2,
         7,
         36,
         38,
         40,
         42,
         44,
         24,
         44,
         12,
         32,
         44,
         32,
         36,
         42,
         34,
         44,
         28],
)

function temperature_impl(mt::SinglePH,model::IF97Region{:r2a},p,_h)
    P = convert_unit(u"Pa",u"MPa",p)
    h = convert_unit(u"J/kg",u"kJ/kg",_h)
    hstar = 2000.0
    Pstar = 1.0
    n = Region2a_t_ph.n
    I = Region2a_t_ph.I 
    J = Region2a_t_ph.J
    π = P / Pstar
    ηη = h / hstar - 2.1

    return sum(n[i]*π^I[i]*ηη^J[i] for i=1:34)
end



const Region2b_t_ph = (
    n = [0.148_950_410_795_16E4,
    0.743_077_983_140_34E3,
    -0.977_083_187_978_37E2,
    0.247_424_647_056_74E1,
    -0.632_813_200_160_26,
    0.113_859_521_296_58E1,
    -0.478_118_636_486_25,
    0.852_081_234_315_44E-2,
    0.937_471_473_779_32,
    0.335_931_186_049_16E1,
    0.338_093_556_014_54E1,
    0.168_445_396_719_04,
    0.738_757_452_366_95,
    -0.471_287_374_361_86,
    0.150_202_731_397_07,
    -0.217_641_142_197_50E-2,
    -0.218_107_553_247_61E-1,
    -0.108_297_844_036_77,
    -0.463_333_246_358_12E-1,
    0.712_803_519_595_51E-4,
    0.110_328_317_899_99E-3,
    0.189_552_483_879_02E-3,
    0.308_915_411_605_37E-2,
    0.135_555_045_549_49E-2,
    0.286_402_374_774_56E-6,
    -0.107_798_573_575_12E-4,
    -0.764_627_124_548_14E-4,
    0.140_523_928_183_16E-4,
    -0.310_838_143_314_34E-4,
    -0.103_027_382_121_03E-5,
    0.282_172_816_350_40E-6,
    0.127_049_022_719_45E-5,
    0.738_033_534_682_92E-7,
    -0.110_301_392_389_09E-7,
    -0.814_563_652_078_33E-13,
    -0.251_805_456_829_62E-10,
    -0.175_652_339_694_07E-17,
    0.869_341_563_441_63E-14],

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
    1,
    1,
    2,
    2,
    2,
    2,
    3,
    3,
    3,
    3,
    4,
    4,
    4,
    4,
    4,
    4,
    5,
    5,
    5,
    6,
    7,
    7,
    9,
    9],

    J = [0,
    1,
    2,
    12,
    18,
    24,
    28,
    40,
    0,
    2,
    6,
    12,
    18,
    24,
    28,
    40,
    2,
    8,
    18,
    40,
    1,
    2,
    12,
    24,
    2,
    12,
    18,
    24,
    28,
    40,
    18,
    24,
    40,
    28,
    2,
    28,
    1,
    40]
)

function temperature_impl(mt::SinglePH,model::IF97Region{:r2b},p,_h)
    P = convert_unit(u"Pa",u"MPa",p)
    h = convert_unit(u"J/kg",u"kJ/kg",_h)
    hstar = 2000.0
    Pstar = 1.0
    n = Region2b_t_ph.n
    I = Region2b_t_ph.I 
    J = Region2b_t_ph.J
    π = P / Pstar
    ηη = h / hstar - 2.1

    return sum(n[i]*π^I[i]*ηη^J[i] for i=1:38)
end

const Region2c_t_ph=(n = [-0.323_683_985_552_42E13,
    0.732_633_509_021_81E13,
    0.358_250_899_454_47E12,
    -0.583_401_318_515_90E12,
    -0.107_830_682_174_70E11,
    0.208_255_445_631_71E11,
    0.610_747_835_645_16E6,
    0.859_777_225_355_80E6,
    -0.257_457_236_041_70E5,
    0.310_810_884_227_14E5,
    0.120_823_158_659_36E4,
    0.482_197_551_092_55E3,
    0.379_660_012_724_86E1,
    -0.108_429_848_800_77E2,
    -0.453_641_726_766_60E-1,
    0.145_591_156_586_98E-12,
    0.112_615_974_072_30E-11,
    -0.178_049_822_406_86E-10,
    0.123_245_796_908_32E-6,
    -0.116_069_211_309_84E-5,
    0.278_463_670_885_54E-4,
    -0.592_700_384_741_76E-3,
    0.129_185_829_918_78E-2],

    I = [-7,
        -7,
        -6,
        -6,
        -5,
        -5,
        -2,
        -2,
        -1,
        -1,
        0,
        0,
        1,
        1,
        2,
        6,
        6,
        6,
        6,
        6,
        6,
        6,
        6],

    J = [0,
        4,
        0,
        2,
        0,
        2,
        0,
        1,
        0,
        2,
        0,
        1,
        4,
        8,
        4,
        0,
        1,
        4,
        10,
        12,
        16,
        20,
        22]
)


function temperature_impl(mt::SinglePH,model::IF97Region{:r2c},p,_h)
    P = convert_unit(u"Pa",u"MPa",p)
    h = convert_unit(u"J/kg",u"kJ/kg",_h)
    hstar = 2000.0
    Pstar = 1.0
    n = Region2c_t_ph.n
    I = Region2c_t_ph.I 
    J = Region2c_t_ph.J
    ππ = P / Pstar + 25
    ηη = h / hstar - 1.8

    return sum(n[i]*ππ^I[i]*ηη^J[i] for i=1:23)
end

const Region2a_t_ps = (n = [-0.392_359_838_619_84E6,
    0.515_265_738_272_70E6,
    0.404_824_431_610_48E5,
    -0.321_937_909_239_02E3,
    0.969_614_242_186_94E2,
    -0.228_678_463_717_73E2,
    -0.449_429_141_243_57E6,
    -0.501_183_360_201_66E4,
    0.356_844_635_600_15,
    0.442_353_358_481_90E5,
    -0.136_733_888_117_08E5,
    0.421_632_602_078_64E6,
    0.225_169_258_374_75E5,
    0.474_421_448_656_46E3,
    -0.149_311_307_976_47E3,
    -0.197_811_263_204_52E6,
    -0.235_543_994_707_60E5,
    -0.190_706_163_020_76E5,
    0.553_756_698_831_64E5,
    0.382_936_914_373_63E4,
    -0.603_918_605_805_67E3,
    0.193_631_026_203_31E4,
    0.426_606_436_986_10E4,
    -0.597_806_388_727_18E4,
    -0.704_014_639_268_62E3,
    0.338_367_841_075_53E3,
    0.208_627_866_351_87E2,
    0.338_341_726_561_96E-1,
    -0.431_244_284_148_93E-4,
    0.166_537_913_564_12E3,
    -0.139_862_920_558_98E3,
    -0.788_495_479_998_72,
    0.721_324_117_538_72E-1,
    -0.597_548_393_982_83E-2,
    -0.121_413_589_539_04E-4,
    0.232_270_967_338_71E-6,
    -0.105_384_635_661_94E2,
    0.207_189_254_965_02E1,
    -0.721_931_552_604_27E-1,
    0.207_498_870_811_20E-6,
    -0.183_406_579_113_79E-1,
    0.290_362_723_486_96E-6,
    0.210_375_278_936_19,
    0.256_812_397_299_99E-3,
    -0.127_990_029_337_81E-1,
    -0.821_981_026_520_18E-5],

    I = [-1.5,
    -1.5,
    -1.5,
    -1.5,
    -1.5,
    -1.5,
    -1.25,
    -1.25,
    -1.25,
    -1.0,
    -1.0,
    -1.0,
    -1.0,
    -1.0,
    -1.0,
    -0.75,
    -0.75,
    -0.5,
    -0.5,
    -0.5,
    -0.5,
    -0.25,
    -0.25,
    -0.25,
    -0.25,
    0.25,
    0.25,
    0.25,
    0.25,
    0.5,
    0.5,
    0.5,
    0.5,
    0.5,
    0.5,
    0.5,
    0.75,
    0.75,
    0.75,
    0.75,
    1.0,
    1.0,
    1.25,
    1.25,
    1.5,
    1.5],

    J = [-24,
    -23,
    -19,
    -13,
    -11,
    -10,
    -19,
    -15,
    -6,
    -26,
    -21,
    -17,
    -16,
    -9,
    -8,
    -15,
    -14,
    -26,
    -13,
    -9,
    -7,
    -27,
    -25,
    -11,
    -6,
    1,
    4,
    8,
    11,
    0,
    1,
    5,
    6,
    10,
    14,
    16,
    0,
    4,
    9,
    17,
    7,
    18,
    3,
    15,
    5,
    18]
)

function temperature_impl(mt::SinglePS,model::IF97Region{:r2a},p,_s)
    P = convert_unit(u"Pa",u"MPa",p)
    s = convert_unit(u"J/(kg*K)",u"kJ/(kg*K)",_s)
    hstar = 2000.0
    Pstar = 1.0
    n = Region2a_t_ps.n
    I = Region2a_t_ps.I 
    J = Region2a_t_ps.J
    π = P / Pstar
    σσ = s / sstar - 2.0

    return sum([n[i]*π^I[i]*σσ^J[i] for i=1:46])
end

const Region2b_t_ps = ( n = [0.316_876_650_834_97E6,
    0.208_641_758_818_58E2,
    -0.398_593_998_035_99E6,
    -0.218_160_585_188_77E2,
    0.223_697_851_942_42E6,
    -0.278_417_034_458_17E4,
    0.992_074_360_714_80E1,
    -0.751_975_122_991_57E5,
    0.297_086_059_511_58E4,
    -0.344_068_785_485_26E1,
    0.388_155_642_491_15,
    0.175_112_950_857_50E5,
    -0.142_371_128_544_49E4,
    0.109_438_033_641_67E1,
    0.899_716_193_084_95,
    -0.337_597_400_989_58E4,
    0.471_628_858_183_55E3,
    -0.191_882_419_936_79E1,
    0.410_785_804_921_96,
    -0.334_653_781_720_97,
    0.138_700_347_775_05E4,
    -0.406_633_261_958_38E3,
    0.417_273_471_596_10E2,
    0.219_325_494_345_32E1,
    -0.103_200_500_090_77E1,
    0.358_829_435_167_03,
    0.525_114_537_260_66E-2,
    0.128_389_164_507_05E2,
    -0.286_424_372_193_81E1,
    0.569_126_836_648_55,
    -0.999_629_545_849_31E-1,
    -0.326_320_377_784_59E-2,
    0.233_209_225_767_23E-3,
    -0.153_348_098_574_50,
    0.290_722_882_399_02E-1,
    0.375_347_027_411_67E-3,
    0.172_966_917_024_11E-2,
    -0.385_560_508_445_04E-3,
    -0.350_177_122_926_08E-4,
    -0.145_663_936_314_92E-4,
    0.564_208_572_672_69E-5,
    0.412_861_500_746_05E-7,
    -0.206_846_711_188_24E-7,
    0.164_093_936_747_25E-8],

    I = [-6,
    -6,
    -5,
    -5,
    -4,
    -4,
    -4,
    -3,
    -3,
    -3,
    -3,
    -2,
    -2,
    -2,
    -2,
    -1,
    -1,
    -1,
    -1,
    -1,
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
    3,
    3,
    3,
    4,
    4,
    5,
    5,
    5],

    J = [0,
    11,
    0,
    11,
    0,
    1,
    11,
    0,
    1,
    11,
    12,
    0,
    1,
    6,
    10,
    0,
    1,
    5,
    8,
    9,
    0,
    1,
    2,
    4,
    5,
    6,
    9,
    0,
    1,
    2,
    3,
    7,
    8,
    0,
    1,
    5,
    0,
    1,
    3,
    0,
    1,
    0,
    1,
    2]
)

function temperature_impl(mt::SinglePS,model::IF97Region{:r2b},p,_s)
    P = convert_unit(u"Pa",u"MPa",p)
    s = convert_unit(u"J/(kg*K)",u"kJ/(kg*K)",_s)
    sstar = 0.7853
    Pstar = 1.0
    n = Region2b_t_ps.n
    I = Region2b_t_ps.I 
    J = Region2b_t_ps.J
    π = P / Pstar
    σσ = 10.0 - s / sstar
    return sum([n[i]*π^I[i]*σσ^J[i] for i=1:44])
end

const Region2c_t_ps = (n = [0.909_685_010_053_65E3,
    0.240_456_670_884_20E4,
    -0.591_623_263_871_30E3,
    0.541_454_041_280_74E3,
    -0.270_983_084_111_92E3,
    0.979_765_250_979_26E3,
    -0.469_667_729_594_35E3,
    0.143_992_746_047_23E2,
    -0.191_042_042_304_29E2,
    0.532_991_671_119_71E1,
    -0.212_529_753_759_34E2,
    -0.311_473_344_137_60,
    0.603_348_408_946_23,
    -0.427_648_397_025_09E-1,
    0.581_855_972_552_59E-2,
    -0.145_970_082_847_53E-1,
    0.566_311_756_310_27E-2,
    -0.761_558_645_845_77E-4,
    0.224_403_429_193_32E-3,
    -0.125_610_950_134_13E-4,
    0.633_231_326_609_34E-6,
    -0.205_419_896_753_75E-5,
    0.364_053_703_900_82E-7,
    -0.297_598_977_892_15E-8,
    0.101_366_185_297_63E-7,
    0.599_257_196_923_51E-11,
    -0.206_778_701_051_64E-10,
    -0.208_742_781_818_86E-10,
    0.101_621_668_250_89E-9,
    -0.164_298_282_813_47E-9],

    I = [-2,
    -2,
    -1,
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
    3,
    3,
    3,
    4,
    4,
    4,
    5,
    5,
    5,
    6,
    6,
    7,
    7,
    7,
    7,
    7],

    J = [0,
    1,
    0,
    0,
    1,
    2,
    3,
    0,
    1,
    3,
    4,
    0,
    1,
    2,
    0,
    1,
    5,
    0,
    1,
    4,
    0,
    1,
    2,
    0,
    1,
    0,
    1,
    3,
    4,
    5]
)

function temperature_impl(mt::SinglePS,model::IF97Region{:r2c},p,_s)
    P = convert_unit(u"Pa",u"MPa",p)
    s = convert_unit(u"J/(kg*K)",u"kJ/(kg*K)",_s)
    sstar = 2.9251
    Pstar = 1.0
    n = Region2c_t_ps.n
    I = Region2c_t_ps.I 
    J = Region2c_t_ps.J
        π = P / Pstar
    σσ = 2.0 - s / sstar
    return sum(n[i]*π^I[i]*σσ^J[i] for i=1:30)
end

