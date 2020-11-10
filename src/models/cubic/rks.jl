struct RedlichKwongSoave{S,T,A} <: CubicModel
    type::S
    tc::T
    pc::T
    ω::T
    mw::T
    vc::Union{T,Nothing}
    _a::T
    _b::T
    aij::A

    function RedlichKwongSoave(tc,pc,ω,mw,vc=nothing,aij = nothing)
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

volume_solver_type(model::RedlichKwongSoave) = CubicRoots()
single_sat_aprox(model::RedlichKwongSoave{SINGLE}) =RKSatAprox(model)
mol_density(model::RedlichKwongSoave{SINGLE},::CriticalPoint,unit=u"mol/(m^3)") = convert_unit(u"mol/L",unit,inv(only(model.vc)))
pressure(model::RedlichKwongSoave{SINGLE},::CriticalPoint,unit=u"Pa") = convert_unit(u"Pa",unit,only(model.pc))
temperature(model::RedlichKwongSoave{SINGLE},::CriticalPoint,unit=u"K") = convert_unit(u"K",unit,only(model.tc))
mol_volume(model::RedlichKwongSoave{SINGLE},::CriticalPoint,unit=u"m^3/mol") = convert_unit(u"m^3/mol",unit,only(model.vc))
acentric_factor(model::RedlichKwongSoave{SINGLE}) = only(model.ω)
molecular_weight(model::RedlichKwongSoave{SINGLE}) = only(model.mw)
molecular_weight(model::RedlichKwongSoave{MULTI}) = model.mw

function mol_density(model::RedlichKwongSoave{MULTI},::CriticalPoint,unit=u"mol/(m^3)")
    return convert_unit.(u"mol/L",unit,1 ./ model.vc)
end

function pressure(model::RedlichKwongSoave{MULTI},::CriticalPoint,unit=u"Pa")
    return convert_unit.(u"Pa",unit,model.pc)
end

function temperature(model::RedlichKwongSoave{MULTI},::CriticalPoint,unit=u"K")
    return convert_unit.(u"K",unit,model.tc)
end

function mol_volume(model::RedlichKwongSoave{MULTI},::CriticalPoint,unit=u"m^3/mol")
    return convert_unit.(u"m^3/mol",unit,model.vc)
end

function acentric_factor(model::RedlichKwongSoave{MULTI})
    return model.ω
end

function RedlichKwongSoave(;tc,pc,ω,mw,vc=nothing,aij=nothing)
    return RedlichKwongSoave(tc,pc,ω,mw,vc,aij)
end

function RedlichKwongSoave(model::ThermoModel)
    tc = temperature(model,CriticalPoint())
    pc = pressure(model,CriticalPoint())
    ω = acentric_factor(model)
    mw = molecular_weight(model)
    return RedlichKwongSoave(tc,pc,ω,mw)
end

function cubic_aα(model::RedlichKwongSoave,t,i,j)
    m_poly = (0.47979, 1.576, 0.1925, 0.025)
    _1 = one(t)
    aᵢ =model._a[i]
    tcᵢ = model.tc[i]
    
    
    if tcᵢ > t #α(t) is 1 for supercritical values
        sqrt_trᵢ = sqrt(t/tcᵢ) 
        ωᵢ = model.ω[i]
        mᵢ = evalpoly(ωᵢ,m_poly)
        sqrt_αᵢ = (_1+mᵢ * √(_1-sqrt_trᵢ))
    else
        sqrt_αᵢ = _1
    end
    if i === j
        return  aᵢ* sqrt_αᵢ^2
    else
        aⱼ =model._a[j]
        tcⱼ = model.tc[j]
        
        if tcⱼ > t #α(t) is 1 for supercritical values
            sqrt_trⱼ = sqrt(t/tcⱼ) 
            ωⱼ = model.ω[j]
            mⱼ = evalpoly(ωⱼ,m_poly)
            sqrt_αⱼ = (_1+mⱼ * √(_1-sqrt_trⱼ))
        else
            sqrt_αⱼ = _1 
        end

        return sqrt(aᵢ*aⱼ)*sqrt_αᵢ*sqrt_αⱼ
    end
end

function cubic_ab(mt::SinglePT,model::RedlichKwongSoave{SINGLE},v,t)
    a = cubic_aα(model,t,1,1)
    b = only(model._b)
    return a,b
end




function cubic_ab(mt::MultiPT,model::RedlichKwongSoave{MULTI},p,t,x)
    #two options to introduce alpha:
    #here: it will allocate, but less ops
    bi = model._b
    b = dot(bi,x)
    #here: it will not allocate, but more ops
    sss= (i,j)->cubic_aα(model,t,i,j)
    a = cubic_mixing_rule(sss, x, model.aij)
    return a,b
end

function cubic_abp(mt::SingleVT,model::RedlichKwongSoave{SINGLE},v,t)
    a,b = cubic_ab(QuickStates.pt(),model,v,t) #v is ignored
    p = RGAS*t/(v-b) - a/(v*v)
    return a,b,p
end

function cubic_abp(mt::MultiVT,model::RedlichKwongSoave{MULTI},v,t,x)
    a,b = cubic_ab(QuickStates.ptx(),model,v,t,x) #v is ignored
    p =  RGAS*t/(v-b) - a/((v+b)*v)
    return a,b,p
end

function fugacity_coeff_impl(mt::SingleVT,model::RedlichKwongSoave{SINGLE},v,t)
    a,b,p =  cubic_abp(mt,model,v,t)
    RTinv = 1/(RGAS*t)
    A = a*p*RTinv*RTinv
    B = b*p*RTinv
    z = p*v*RTinv
     _1 = one(z)
     logϕ = z - _1 - log(z-B) - A/z
end

function  αR_impl(mt::MultiVT,model::RedlichKwongSoave{MULTI},rho,t,x)
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
function cubic_poly(mt::SinglePT,model::RedlichKwongSoave{SINGLE},p,t)
    a,b = cubic_ab(QuickStates.pt(),model,p,t)
    RTinv = 1/(RGAS*t)
    A = a*p*RTinv*RTinv
    B = b*p*RTinv
    _1 = one(a)
    return (-A*B, -B*(B+_1) + A, -_1, _1)
end

function cubic_poly(mt::MultiPT,model::RedlichKwongSoave{MULTI},p,t,x)
    a,b = cubic_ab(QuickStates.ptx(),model,p,t,x)
    RTinv = 1/(RGAS*t)
    A = a*p*RTinv*RTinv
    B = b*p*RTinv
    _1 = one(a)
    return (-A*B, -B*(B+_1) + A, -_1, _1)
end


