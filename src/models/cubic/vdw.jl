struct VanDerWaals{S,T,A} <: CubicModel
    type::S
    tc::T
    pc::T
    ω::T
    mw::T
    vc::Union{T,Nothing}
    _a::T
    _b::T
    aij::A

    function VanDerWaals(tc,pc,ω,mw,vc=nothing,aij = nothing)
        if length(tc) == 1
            type = SINGLE()

        else 
            type = MULTI()
        end
        S = typeof(type)
        tc,pc,ω,mw = promote(tc,pc,ω,mw)
        T = typeof(tc)
        A = typeof(aij)
        _a = 0.421875 .*((RGAS .*tc).^2) ./pc
        _b = 0.125 .*((RGAS .*tc)) ./pc
        return new{S,T,A}(type,tc,pc,ω,mw,vc,_a,_b,aij)
    end
end

volume_solver_type(model::VanDerWaals) = CubicRoots()
single_sat_Approx(model::VanDerWaals{SINGLE}) =VdWSatApprox(model)
mol_density(model::VanDerWaals{SINGLE},::CriticalPoint,unit=u"mol/(m^3)") = convert_unit(u"mol/L",unit,inv(only(model.vc)))
pressure(model::VanDerWaals{SINGLE},::CriticalPoint,unit=u"Pa") = convert_unit(u"Pa",unit,only(model.pc))
temperature(model::VanDerWaals{SINGLE},::CriticalPoint,unit=u"K") = convert_unit(u"K",unit,only(model.tc))
mol_volume(model::VanDerWaals{SINGLE},::CriticalPoint,unit=u"m^3/mol") = convert_unit(u"m^3/mol",unit,only(model.vc))
acentric_factor(model::VanDerWaals{SINGLE}) = only(model.ω)
molecular_weight(model::VanDerWaals{SINGLE}) = only(model.mw)
molecular_weight(model::VanDerWaals{MULTI}) = model.mw


function mol_density(model::VanDerWaals{MULTI},::CriticalPoint,unit=u"mol/(m^3)")
    return convert_unit.(u"mol/L",unit,1 ./ model.vc)
end

function pressure(model::VanDerWaals{MULTI},::CriticalPoint,unit=u"Pa")
    return convert_unit.(u"Pa",unit,model.pc)
end

function temperature(model::VanDerWaals{MULTI},::CriticalPoint,unit=u"K")
    return convert_unit.(u"K",unit,model.tc)
end

function mol_volume(model::VanDerWaals{MULTI},::CriticalPoint,unit=u"m^3/mol")
    return convert_unit.(u"m^3/mol",unit,model.vc)
end

function acentric_factor(model::VanDerWaals{MULTI})
    return model.ω
end

function VanDerWaals(;tc,pc,ω,mw,vc=nothing,aij=nothing)
    return VanDerWaals(tc,pc,ω,mw,vc,aij)
end

function VanDerWaals(model::ThermoModel)
    tc = temperature(model,CriticalPoint())
    pc = pressure(model,CriticalPoint())
    ω = acentric_factor(model)
    mw = molecular_weight(model)
    return VanDerWaals(tc,pc,ω,mw)
end


function cubic_ab(mt::SinglePT,model::VanDerWaals{SINGLE},v,t)
    a = only(model._a)
    b = only(model._b)
    return a,b
end

function cubic_ab(mt::MultiPT,model::VanDerWaals{MULTI},p,t,x)
    
    #two options to introduce alpha:
    #here: it will allocate, but less ops
    ai = model._a
    bi = model._b
    b = dot(bi,x)
    #here: it will not allocate, but more ops
    mixrule= (ai,aj)->sqrt(ai*aj)
    a = mixing_rule(mixrule, x, ai,model.aij)
    return a,b
end

function cubic_abp(mt::SingleVT,model::VanDerWaals{SINGLE},v,t)
    a,b = cubic_ab(QuickStates.pt(),model,v,t) #v is ignored
    p = RGAS*t/(v-b) - a/(v*v)
    return a,b,p
end

function cubic_abp(mt::MultiVT,model::VanDerWaals{MULTI},v,t,x)
    a,b = cubic_ab(QuickStates.ptx(),model,v,t,x) #v is ignored
    p =  RGAS*t/(v-b) - a/(v*v)
    return a,b,p
end

function fugacity_coeff_impl(mt::SingleVT,model::VanDerWaals{SINGLE},v,t)
    a,b,p =  cubic_abp(mt,model,v,t)
    RTinv = 1/(RGAS*t)
    A = a*p*RTinv*RTinv
    B = b*p*RTinv
    z = p*v*RTinv
     _1 = one(z)
     logϕ = z - _1 - log(z-B) - A/z
end

function  αR_impl(mt::MultiVT,model::VanDerWaals{MULTI},rho,t,x)
    R = RGAS
    RTinv = 1/(RGAS*t)
    v = inv(rho)
    a,b,p =  cubic_abp(mt,model,v,t,x)
    -log(1-b*rho) - a*rho*RTinv
end


function cubic_poly(mt::SinglePT,model::VanDerWaals{SINGLE},p,t)
    a,b = cubic_ab(QuickStates.pt(),model,p,t)
    RTinv = 1/(RGAS*t)
    A = a*p*RTinv*RTinv
    B = b*p*RTinv
    _1 = one(a)
    return (-A*B, A, -B-_1, _1)
end

function cubic_poly(mt::MultiPT,model::VanDerWaals{MULTI},p,t,x)
    a,b = cubic_ab(QuickStates.ptx(),model,p,t,x)
    RTinv = 1/(RGAS*t)
    A = a*p*RTinv*RTinv
    B = b*p*RTinv
    _1 = one(a)
    return (-A*B, A, -B-_1, _1)
end



struct VdWSatApprox{M} <: SaturationModel
    model::M
    function VdWSatApprox(model::CubicModel)
        T = typeof(model)
        return new{T}(model)
    end
end

function VdWSatApprox(tc,pc,mw=copy(pc),vc=nothing,ω=nothing,aij = nothing)
    model = VanDerWaals(tc,pc,mw,vc,ω,aij)
    return VdWSatApprox(model)
end

function VdWSatApprox(;tc,pc,mw=copy(pc),vc=nothing,ω=nothing,aij=nothing)
    model =  VanDerWaals(tc,pc,mw,vc,ω,aij)
    return VdWSatApprox(model)
end

function pressure_impl(mt::SingleSatT,model::VdWSatApprox,t)
    a = 4.406664258927600 
    b = 2.205610041969020 
    c = 0.243757663628277 
    d = -0.206185849671953 
    e = -4.650419726136550 
    f = 0.019582516700758
    tc = only(model.model.tc)
    pc = only(model.model.pc)
    tr = t/tc
    atr = tc/t
    num = evalpoly(atr,(a,c,e)) 
    denom = evalpoly(atr,(1.0,b,d,f))
    pr = exp(num/denom)*tr
    p = pr*pc
end

initial_temperature(model::VdWSatApprox,p) = critical_sat_interpolation(model,p)
