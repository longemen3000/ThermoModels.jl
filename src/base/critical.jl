struct CriticalData{T} <: ThermoModel
    tc::T
    vc::T
    pc::T
    zc::T
    ω::T
end
"""
    CriticalData(;tc = NaN,vc=NaN,pc=NaN,zc=NaN,ω=NaN)

Creates a CriticalData struct, storing information about a critical point.

all default units are/are converted to standard SI units.

If only 3 values are given, the fourth one will be calculated via Pv = ZRT

"""
function CriticalData(;tc = NaN,vc=NaN,pc=NaN,zc=NaN,ω=NaN)
    #try to calculate remaining values
    if !isnan(tc) & !isnan(vc) & !isnan(pc) & isnan(zc)
        _tc = ThermoState.normalize_units(tc)
        _vc = ThermoState.normalize_units(vc)
        _pc = ThermoState.normalize_units(pc)
        _zc  =_pc*_vc/(RGAS*_tc)
        
    elseif isnan(tc) & !isnan(vc) & !isnan(pc) & !isnan(zc)
        _zc = zc
        _vc = ThermoState.normalize_units(vc)
        _pc = ThermoState.normalize_units(pc)
        _tc  =_pc*_vc/(RGAS*_zc)
    elseif !isnan(tc) & isnan(vc) & !isnan(pc) & !isnan(zc)
        _tc = ThermoState.normalize_units(tc)
        _zc = zc
        _pc = ThermoState.normalize_units(pc)
        _vc  =_zc*RGAS*_tc/_pc
    elseif !isnan(tc) & !isnan(vc) & isnan(pc) & !isnan(zc)
        _tc = ThermoState.normalize_units(tc)
        _vc = ThermoState.normalize_units(vc)
        _zc = zc
        _pc  =_zc*RGAS*_tc/_vc
    else
        _tc = ThermoState.normalize_units(tc)
        _vc = ThermoState.normalize_units(vc)
        _pc = ThermoState.normalize_units(pc)
        _zc  =zc
    end
    return CriticalData(_tc,_vc,_pc,_zc,ω)
end

function temperature(model::CriticalData,st::CriticalPoint,unit=u"K")
    return convert_unit(u"K",unit,model.tc)
end

function pressure(model::CriticalData,st::CriticalPoint,unit=u"Pa")
    return convert_unit(u"Pa",unit,model.pc)
end

function mol_volume(model::CriticalData,st::CriticalPoint,unit=u"m^3/mol",mw=nothing)
    return convert_unit(u"m^3/mol",unit,model.vc)
end

function mol_density(model::CriticalData,st::CriticalPoint,unit=u"mol/(m^3)",mw=nothing)
    return convert_unit(u"mol/(m^3)",unit,inv(model.vc))
end

function mass_volume(model::CriticalData,st::CriticalPoint,unit=u"m^3/kg",mw=nothing)
    res = ThermoState.mw_div(model.vc,mw)
    return convert_unit(u"m^3/kg",unit,res)
end

function mass_density(model::CriticalData,st::CriticalPoint,unit=u"kg/(m^3)",mw=nothing)
    res = inv(ThermoState.mw_div(model.vc,mw))
    return convert_unit(u"kg/(m^3)",unit,res)
end

