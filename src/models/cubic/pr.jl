struct PengRobinson{S,T,A} <: CubicModel
    type::S
    tc::T
    pc::T
    ω::T
    mw::T
    vc::Union{T,Nothing}
    _a::T
    _b::T
    aij::A

    function PengRobinson(tc,pc,ω,mw,vc=nothing,aij = nothing)
        if length(tc) == 1
            type = SINGLE()
        else 
            type = MULTI()
        end
        S = typeof(type)
        tc,pc,ω,mw = promote(tc,pc,ω,mw)
        T = typeof(tc)
        A = typeof(aij)
        _a = 0.45723553 .*((RGAS .*tc).^2) ./pc
        _b = 0.07779607 .*((RGAS .*tc)) ./pc
        return new{S,T,A}(type,tc,pc,ω,mw,vc,_a,_b,aij)
    end
end

volume_solver_type(model::PengRobinson) = CubicRoots()
single_sat_Approx(model::PengRobinson{SINGLE}) =PRSatApprox(model)
mol_density(model::PengRobinson{SINGLE},::CriticalPoint,unit=u"mol/(m^3)") = convert_unit(u"mol/L",unit,inv(only(model.vc)))
pressure(model::PengRobinson{SINGLE},::CriticalPoint,unit=u"Pa") = convert_unit(u"Pa",unit,only(model.pc))
temperature(model::PengRobinson{SINGLE},::CriticalPoint,unit=u"K") = convert_unit(u"K",unit,only(model.tc))
mol_volume(model::PengRobinson{SINGLE},::CriticalPoint,unit=u"m^3/mol") = convert_unit(u"m^3/mol",unit,only(model.vc))
acentric_factor(model::PengRobinson{SINGLE}) = only(model.ω)
molecular_weight(model::PengRobinson{SINGLE}) = only(model.mw)
molecular_weight(model::PengRobinson{MULTI}) = model.mw

function mol_density(model::PengRobinson{MULTI},::CriticalPoint,unit=u"mol/(m^3)")
    return convert_unit.(u"mol/L",unit,1 ./ model.vc)
end

function pressure(model::PengRobinson{MULTI},::CriticalPoint,unit=u"Pa")
    return convert_unit.(u"Pa",unit,model.pc)
end

function temperature(model::PengRobinson{MULTI},::CriticalPoint,unit=u"K")
    return convert_unit.(u"K",unit,model.tc)
end

function mol_volume(model::PengRobinson{MULTI},::CriticalPoint,unit=u"m^3/mol")
    return convert_unit.(u"m^3/mol",unit,model.vc)
end

function acentric_factor(model::PengRobinson{MULTI})
    return model.ω
end

function PengRobinson(;tc,pc,ω,mw,vc=nothing,aij=nothing)
    return PengRobinson(tc,pc,ω,mw,vc,aij)
end

function PengRobinson(model::ThermoModel)
    tc = temperature(model,CriticalPoint())
    pc = pressure(model,CriticalPoint())
    ω = acentric_factor(model)
    mw = molecular_weight(model)
    return PengRobinson(tc,pc,ω,mw)
end

function cubic_aα(model::PengRobinson,t,i,j)
    
    m_poly = (0.37464,1.54226,0.26992)
    _1 = one(t)
    aᵢ =model._a[i]
    tcᵢ = model.tc[i]
    sqrt_trᵢ = min(sqrt(t/tcᵢ),_1) #α(t) is one for supercritical values
    ωᵢ = model.ω[i]
    mᵢ = evalpoly(ωᵢ,m_poly)
    if i === j
        return  aᵢ* (_1+mᵢ * √(_1-sqrt_trᵢ))^2
    else
        aⱼ =model._a[j]
        tcⱼ = model.tc[j]
        sqrt_trⱼ = min(sqrt(t/tcⱼ),_1)
        ωⱼ = model.ω[j]
        mⱼ = evalpoly(ωⱼ,m_poly)
        sqrt_αᵢ= (_1+mᵢ * √(_1-sqrt_trᵢ))
        sqrt_αⱼ= (_1+mⱼ * √(_1-sqrt_trⱼ))
        return sqrt(aᵢ*aⱼ)*sqrt_αᵢ*sqrt_αⱼ
    end
end


function cubic_ab(mt::SinglePT,model::PengRobinson{SINGLE},v,t)
    a = cubic_aα(model,t,1,1)
    b = only(model._b)
    return a,b
end

function cubic_ab(mt::MultiPT,model::PengRobinson{MULTI},p,t,x)
    #two options to introduce alpha:
    #here: it will allocate, but less ops
    bi = model._b
    b = dot(bi,x)
    #here: it will not allocate, but more ops
    sss= (i,j)->cubic_aα(model,t,i,j)
    a = cubic_mixing_rule(sss, x, model.aij)
    return a,b
end

function cubic_abp(mt::SingleVT,model::PengRobinson{SINGLE},v,t)
    a,b = cubic_ab(QuickStates.pt(),model,v,t) #v is ignored
    _1 = one(b)
    denom = evalpoly(v,(-b*b,2*b,_1))
    p = RGAS*t/(v-b) - a/denom   
    return a,b,p
end

function cubic_abp(mt::MultiVT,model::PengRobinson{MULTI},v,t,x)
    a,b = cubic_ab(QuickStates.ptx(),model,v,t,x) #v is ignored
    _1 = one(b)
    denom = evalpoly(v,(-b*b,2*b,_1))
    p = RGAS*t/(v-b) - a/denom   
    return a,b,p
end

function fugacity_coeff_impl(mt::SingleVT,model::PengRobinson{SINGLE},v,t)
    a,b,p =  cubic_abp(mt,model,v,t)
    RTinv = 1/(RGAS*t)
    A = a*p*RTinv*RTinv
    B = b*p*RTinv
    z = p*v*RTinv
     _1 = one(z)
     logϕ = z - _1 - log(z-B) - A/z
end

const PRΔ1 = 1+√2 
const PRΔ2 = 1-√2
const ΔPRΔ = 2*√2

function  αR_impl(mt::MultiVT,model::PengRobinson{MULTI},rho,t,x)
    R = RGAS
    RTinv = 1/(RGAS*t)
    v = inv(rho)
    a,b,p =  cubic_abp(mt,model,v,t,x)
    -log(1-b*rho) - a*RTinv*log((PRΔ1*b*rho+1)/(PRΔ2*b*rho+1))/(ΔPRΔ*b)
end


function cubic_poly(mt::SinglePT,model::PengRobinson{SINGLE},p,t)
    a,b = cubic_ab(QuickStates.pt(),model,p,t)
    RTinv = 1/(RGAS*t)
    A = a*p*RTinv*RTinv
    B = b*p*RTinv
    k0 = B*(B*(B+1.0)-A)
    k1 = -B*(2*B+1.0) + A
    k2 = -1.0
    k3 = 1.0
    return (k0,k1,k2,k3)
end

function cubic_poly(mt::MultiPT,model::PengRobinson{MULTI},p,t,x)
    a,b = cubic_ab(QuickStates.ptx(),model,p,t,x)
    RTinv = 1/(RGAS*t)
    A = a*p*RTinv*RTinv
    B = b*p*RTinv
    k0 = B*(B*(B+1.0)-A)
    k1 = -B*(2*B+1.0) + A
    k2 = -1.0
    k3 = 1.0
    return (k0,k1,k2,k3)
end


struct PRSatApprox{M} <: SaturationModel
    model::M
    function PRSatApprox(model::CubicModel)
        T = typeof(model)
        return new{T}(model)
    end
end

function PRSatApprox(tc,pc,mw=copy(pc),vc=nothing,ω=nothing,aij = nothing)
    model = PengRobinson(tc,pc,mw,vc,ω,aij)
    return PRSatApprox(model)
end

function PRSatApprox(;tc,pc,mw=copy(pc),vc=nothing,ω=nothing,aij=nothing)
    model =  PengRobinson(tc,pc,mw,vc,ω,aij)
    return PRSatApprox(model)
end

function pressure_impl(mt::SingleSatT,model::PRSatApprox,t)
    a = -2.605272488488440
    b = -9.017571450539830
    c = -19.896014683288000
    d = 0.579677284391001
    e = 22.501011849124900
    f = -0.041440165126830
    tc = only(model.model.tc)
    pc = only(model.model.pc)
    tr = t/tc
    atr = tc/t
    num = evalpoly(atr,(a,c,e)) 
    denom = evalpoly(atr,(1.0,b,d,f))
    pr = exp(num/denom)*tr
    p = pr*pc
end

initial_temperature(model::PRSatApprox,p) = critical_sat_interpolation(model,p)
