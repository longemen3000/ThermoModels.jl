
initial_temperature(model::SaturationModel,p) = 400.0


function temperature_impl(mt::SingleSatP,model::T,p) where T<:SaturationModel
    f0(t) = pressure_impl(QuickStates.sat_t(),model,t)-p
    df0(t) = ForwardDiff.derivative(f0,t)
    t0 = initial_temperature(model,p)
    t = Roots.find_zero(f0, t0)
    return t
end

temperature_impl(model::SaturationModel,p) =temperature_impl(QuickStates.sat_p(),model,p)

function temperature_impl(mt::SingleΦP,model::SaturationModel,p)
    f0(t) = temperature_impl(QuickStates.sat_t(),model,t)-p
    df0(t) = ForwardDiff.derivative(f0,t)
    t0 = initial_temperature(model,p)
    t = Roots.find_zero((f0, df0), t0, Roots.Newton())
    return t
end

function pressure_impl(mt::SingleΦT,model::SaturationModel,p)
    return pressure_impl(QuickStates.sat_t(),model,t)
end

pressure_impl(model::SaturationModel,t) =pressure_impl(QuickStates.sat_t(),model,t)

function pressure(model::SaturationModel,st::ThermodynamicState,unit=u"Pa")
    return pressure(state_type(st),model,st,unit)
end
function pressure(mt::SingleSatT,model::SaturationModel,st::ThermodynamicState,unit)
    t = temperature(FromState(),st)
    p = pressure_impl(mt,model,t)
    return convert_unit(u"Pa",unit,p)
end

function temperature(model::SaturationModel,st::ThermodynamicState,unit=u"K")
    return temperature(state_type(st),model,st,unit)
end
function temperature(mt::SingleSatP,model::SaturationModel,st::ThermodynamicState,unit)
    p = pressure(FromState(),st)
    t = temperature_impl(mt,model,p)
    return convert_unit(u"K",unit,t)
end

function critical_sat_interpolation(model::SaturationModel,p)
    #ln(pr) = h*(1-1/tr)
    #1/(1-ln(pr)/h)*tc = t
    tc = temperature(model,CriticalPoint())
    pc = pressure(model,CriticalPoint())
    t7 = 0.7*tc
    p7 = pressure_impl(QuickStates.sat_t(),model,t7)
    
    h = 2.3333333333333335*log(pc/p7)
    return 1/(1-log(p/pc)/h)*tc
end
