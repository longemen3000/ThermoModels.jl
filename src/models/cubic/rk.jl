struct RedlichKwong{S,T,A} <: CubicModel
    type::S
    tc::T
    pc::T
    ω::T
    mw::T
    vc::Union{T,Nothing}
    _a::T
    _b::T
    aij::A

    function RedlichKwong(tc,pc,ω,mw,vc=nothing,aij = nothing)
        if length(tc) == 1
            type = SINGLE()

        else 
            type = MULTI()
        end
        S = typeof(type)
        tc,pc,ω,mw = promote(tc,pc,ω,mw)
        T = typeof(tc)
        A = typeof(aij)
        _a = 0.42748 .*((RGAS .*tc).^2) ./pc
        _b = 0.08664 .*((RGAS .*tc)) ./pc
        return new{S,T,A}(type,tc,pc,ω,mw,vc,_a,_b,aij)
    end
end

volume_solver_type(model::RedlichKwong) = CubicRoots()
single_sat_Approx(model::RedlichKwong{SINGLE}) =RKSatApprox(model)
mol_density(model::RedlichKwong{SINGLE},::CriticalPoint,unit=u"mol/(m^3)") = convert_unit(u"mol/L",unit,inv(only(model.vc)))
pressure(model::RedlichKwong{SINGLE},::CriticalPoint,unit=u"Pa") = convert_unit(u"Pa",unit,only(model.pc))
temperature(model::RedlichKwong{SINGLE},::CriticalPoint,unit=u"K") = convert_unit(u"K",unit,only(model.tc))
mol_volume(model::RedlichKwong{SINGLE},::CriticalPoint,unit=u"m^3/mol") = convert_unit(u"m^3/mol",unit,only(model.vc))
acentric_factor(model::RedlichKwong{SINGLE}) = only(model.ω)
molecular_weight(model::RedlichKwong{SINGLE}) = only(model.mw)
molecular_weight(model::RedlichKwong{MULTI}) = model.mw

function mol_density(model::RedlichKwong{MULTI},::CriticalPoint,unit=u"mol/(m^3)")
    return convert_unit.(u"mol/L",unit,1 ./ model.vc)
end

function pressure(model::RedlichKwong{MULTI},::CriticalPoint,unit=u"Pa")
    return convert_unit.(u"Pa",unit,model.pc)
end

function temperature(model::RedlichKwong{MULTI},::CriticalPoint,unit=u"K")
    return convert_unit.(u"K",unit,model.tc)
end

function mol_volume(model::RedlichKwong{MULTI},::CriticalPoint,unit=u"m^3/mol")
    return convert_unit.(u"m^3/mol",unit,model.vc)
end

function acentric_factor(model::RedlichKwong{MULTI})
    return model.ω
end

function RedlichKwong(;tc,pc,ω,mw,vc=nothing,aij=nothing)
    return RedlichKwong(tc,pc,ω,mw,vc,aij)
end

function RedlichKwong(model::ThermoModel)
    tc = temperature(model,CriticalPoint())
    pc = pressure(model,CriticalPoint())
    ω = acentric_factor(model)
    mw = molecular_weight(model)
    return RedlichKwong(tc,pc,ω,mw)
end

function cubic_aα(model::RedlichKwong,t,i,j)
    aᵢ =model._a[i]
    tcᵢ = model.tc[i]
    sqrt_t⁻¹ = sqrt(inv(t))
    if i === j
        return  aᵢ* sqrt(tcᵢ)*sqrt_t⁻¹
    else
        aⱼ =model._a[j]
        tcⱼ = model.tc[j]
        return sqrt(aᵢ*aⱼ* sqrt(tcᵢ* tcⱼ))*sqrt_t⁻¹
    end
end

function cubic_ab(mt::SinglePT,model::RedlichKwong{SINGLE},v,t)
    a = cubic_aα(model,t,1,1)
    b = only(model._b)
    return a,b
end




function cubic_ab(mt::MultiPT,model::RedlichKwong{MULTI},p,t,x)
    #two options to introduce alpha:
    #here: it will allocate, but less ops
    bi = model._b
    b = dot(bi,x)
    #here: it will not allocate, but more ops
    sss= (i,j)->cubic_aα(model,t,i,j)
    a = cubic_mixing_rule(sss, x, model.aij)
    return a,b
end

function cubic_abp(mt::SingleVT,model::RedlichKwong{SINGLE},v,t)
    a,b = cubic_ab(QuickStates.pt(),model,v,t) #v is ignored
    p = RGAS*t/(v-b) - a/(v*v)
    return a,b,p
end

function cubic_abp(mt::MultiVT,model::RedlichKwong{MULTI},v,t,x)
    a,b = cubic_ab(QuickStates.ptx(),model,v,t,x) #v is ignored
    p =  RGAS*t/(v-b) - a/((v+b)*v)
    return a,b,p
end

function fugacity_coeff_impl(mt::SingleVT,model::RedlichKwong{SINGLE},v,t)
    a,b,p =  cubic_abp(mt,model,v,t)
    RTinv = 1/(RGAS*t)
    A = a*p*RTinv*RTinv
    B = b*p*RTinv
    z = p*v*RTinv
     _1 = one(z)
     logϕ = z - _1 - log(z-B) - A/z
end

function  αR_impl(mt::MultiVT,model::RedlichKwong{MULTI},rho,t,x)
    R = RGAS
    RTinv = 1/(RGAS*t)
    v = inv(rho)
    a,b,p =  cubic_abp(mt,model,v,t,x)
    -log(1-b*rho) - a*RTinv*log(b*rho+1)/b
end
    #=
    k0 = -AB
    k1 = -B*(B+1) + A
    k2 = -1
    k3 = 1
    =#
function cubic_poly(mt::SinglePT,model::RedlichKwong{SINGLE},p,t)
    a,b = cubic_ab(QuickStates.pt(),model,p,t)
    RTinv = 1/(RGAS*t)
    A = a*p*RTinv*RTinv
    B = b*p*RTinv
    _1 = one(a)
    return (-A*B, -B*(B+_1) + A, -_1, _1)
end

function cubic_poly(mt::MultiPT,model::RedlichKwong{MULTI},p,t,x)
    a,b = cubic_ab(QuickStates.ptx(),model,p,t,x)
    RTinv = 1/(RGAS*t)
    A = a*p*RTinv*RTinv
    B = b*p*RTinv
    _1 = one(a)
    return (-A*B, -B*(B+_1) + A, -_1, _1)
end


struct RKSatApprox{M} <: SaturationModel
    model::M
    function RKSatApprox(model::CubicModel)
        T = typeof(model)
        return new{T}(model)
    end
end

function RKSatApprox(tc,pc,mw=copy(pc),vc=nothing,ω=nothing,aij = nothing)
    model = RedlichKwong(tc,pc,mw,vc,ω,aij)
    return RKSatApprox(model)
end

function RKSatApprox(;tc,pc,mw=copy(pc),vc=nothing,ω=nothing,aij=nothing)
    model =  RedlichKwong(tc,pc,mw,vc,ω,aij)
    return RKSatApprox(model)
end

function pressure_impl(mt::SingleSatT,model::RKSatApprox,t)
    a = 4.869351869381000
    b = 2.783631394064670
    c = 1.106187857326840
    d = -0.256386265669017
    e = -5.975540029435890
    f = 0.024270608586173
    tc = only(model.model.tc)
    pc = only(model.model.pc)
    tr = t/tc
    atr = tc/t
    num = evalpoly(atr,(a,c,e)) 
    denom = evalpoly(atr,(1.0,b,d,f))
    pr = exp(num/denom)*tr
    p = pr*pc
end

initial_temperature(model::RKSatApprox,p) = critical_sat_interpolation(model,p)
