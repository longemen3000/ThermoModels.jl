struct VanDerWaals{S,T,A} <: CubicModel
    type::S
    tc::T
    pc::T
    vc::T
    ω::T
    mw::T
    _a::T
    _b::T
    aij::A
end

function VanDerWaals(tc,pc,vc,ω,mw,aij = nothing)
    if length(tc) == 1
        type = SingleComponent()
        _a = 0.421875*((RGAS*tc)^2)/pc
        _b = 0.125*((RGAS*tc))/pc
    else 
        type = MaterialCompounds{MOL,FRACTION}()
        _a = 0.421875 .*((RGAS .*tc).^2) ./pc
        _b = 0.125 .*((RGAS .*tc)) ./pc
    end
    
    return VanDerWaals(type,tc,pc,vc,ω,mw,_a,_b,aij)
end

VanDerWaals(;tc,pc,vc,ω,mw) = VanDerWaals(tc,pc,vc,ω,mw)

function VanDerWaals(model::ThermoModel)
    tc = temperature(model,CriticalPoint())
    pc = pressure(model,CriticalPoint())
    ω = acentric_factor(model)
    vc = mol_volume(model,CriticalPoint())
    mw = molecular_weight(model)
    return VanDerWaals(tc,pc,vc,ω,mw)
end


function pressure_impl(mt::SingleVT,model::VanDerWaals,v,t)
    return RGAS*t/(v-model._b) - model._a/(v*v)
end
