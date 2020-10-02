function mol_volume_impl(mt::SingleΦT,model::SatLiquidModel,t)
    return mol_volume_impl(QuickStates.sat_t(),model,t)
end

function mol_volume_impl(mt::MultiΦT,model::SatLiquidModel,t,x)
    return mol_volume_impl(QuickStates.sat_tx(),model,t,x)
end

function mol_volume(model::SatLiquidModel,st::ThermodynamicState,unit=u"m^3/mol",mw=nothing)
    return mol_volume(state_type(st),model,st,unit)
end


function mass_volume(model::SatLiquidModel,st::ThermodynamicState,unit=u"m^3/kg",mw=nothing)
    return mol_volume(state_type(st),model,st,unit)
end

function total_volume(model::SatLiquidModel,st::ThermodynamicState,unit=u"m^3",mw=nothing)
    return total_volume(state_type(st),model,st,unit)
end

function mol_density(model::SatLiquidModel,st::ThermodynamicState,unit=u"mol/m^3",mw=nothing)
    return mol_density(state_type(st),model,st,unit)
end

function mass_density(model::SatLiquidModel,st::ThermodynamicState,unit=u"kg/m^3",mw=nothing)
    return mass_density(state_type(st),model,st,unit)
end

function mol_volume(mt::SingleSatT,model::SatLiquidModel,st::ThermodynamicState,unit,mw)
    t = temperature(FromState(),st) 
    res = mol_volume_impl(mt,model,t)
    return convert_unit(u"m^3/mol",unit,res)
end

function mol_volume(mt::MultiSatT,model::SatLiquidModel,st::ThermodynamicState,unit,mw)
    t = temperature(FromState(),st) 
    x = mol_fraction(FromState(),st,nothing,mw) 
    res = mol_volume_impl(mt,model,t,x)
    return convert_unit(u"m^3/mol",unit,res)
end

function mass_volume(mt::SingleSatT,model::SatLiquidModel,st::ThermodynamicState,unit,mw)
    t = temperature(FromState(),st) 
    res = mol_volume_impl(mt,model,t)
    return convert_unit(u"m^3/kg",unit,res)
end

function mass_volume(mt::MultiSatT,model::SatLiquidModel,st::ThermodynamicState,unit,mw)
    t = temperature(FromState(),st) 
    x = mol_fraction(FromState(),st,nothing,mw) 
    res = mol_volume_impl(mt,model,t,x)/molar_mass(FromState(),st,u"kg/mol",mw)
    return convert_unit(u"m^3/kg",unit,)
end

function total_volume(mt::SingleSatT,model::SatLiquidModel,st::ThermodynamicState,unit,mw)
    t = temperature(FromState(),st) 
    res = mol_volume_impl(mt,model,t)
    return convert_unit(u"m^3",unit,res)
end

function total_volume(mt::MultiSatT,model::SatLiquidModel,st::ThermodynamicState,unit,mw)
    t = temperature(FromState(),st) 
    x = mol_fraction(FromState(),st,nothing,mw) 
    res = mol_volume_impl(mt,model,t,x)*moles(FromState(),st,u"mol",mw)
    return convert_unit(u"m^3",unit,)
end

function mol_density(mt::SingleSatT,model::SatLiquidModel,st::ThermodynamicState,unit,mw)
    t = temperature(FromState(),st) 
    res = mol_volume_impl(mt,model,t)
    return convert_unit(u"mol/m^3",unit,one(res)/res)
end

function mol_density(mt::MultiSatT,model::SatLiquidModel,st::ThermodynamicState,unit,mw)
    t = temperature(FromState(),st) 
    x = mol_fraction(FromState(),st,nothing,mw) 
    res = mol_volume_impl(mt,model,t,x)
    return convert_unit(u"mol/m^3",unit,one(res)/res)
end

function mass_density(mt::SingleSatT,model::SatLiquidModel,st::ThermodynamicState,unit,mw)
    t = temperature(FromState(),st) 
    res = mol_volume_impl(mt,model,t)
    return convert_unit(u"kg/m^3",unit,one(res)/res)
end

function mass_density(mt::MultiSatT,model::SatLiquidModel,st::ThermodynamicState,unit,mw)
    t = temperature(FromState(),st) 
    x = mol_fraction(FromState(),st,nothing,mw) 
    res = mol_volume_impl(mt,model,t,x)/molar_mass(FromState(),st,u"kg/mol",mw)
    return convert_unit(u"kg/m^3",unit,one(res)/res)
end

