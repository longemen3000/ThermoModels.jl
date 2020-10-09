struct MixRule{M<:ThermoModel,X} <: ThermoModel
    model::M
    xn::X
end

function mol_fraction(st::MixRule{M,X}) where {M,X}
    return st.xn
end

function thermomodel(st::MixRule{M,X}) where {M,X}
    return st.model
end

abstract type AbstractMixtureRule end

"""
    LinarMixing(p)

corresponds to the mix rule Σxᵢpᵢ.

""" 
struct LinearMixing{M} <: AbstractMixtureRule
    prop::M
end

function mixrule(model::T,x) where T<: LinearMixing
    return dot(model.prop,x)
end
"""
    CuadraticMixing(p,f=(x,y)/2,Aij= nothing)

corresponds to the mix rule ΣAᵢⱼxᵢxⱼf(pᵢ,pⱼ)
where f is an averaging function. That means:
- f(pᵢ,pⱼ) = f(pⱼ,pᵢ)
- f(pₖ,pₖ) = pₖ

`Aij` is the weighed interaction matrix. if not provided, A[i,j] == Aᵢⱼ == 1 
""" 
struct QuadraticMixing{P,F,M} <: AbstractMixtureRule
    prop::P
    func::F
    Aij::M
end


function QuadraticMixing(prop)
    return QuadraticMixing(prop,(x,y)->0.5*(x+y),nothing)
end


function QuadraticMixing(prop,f)
    return QuadraticMixing(prop,f,nothing)
end


function mixrule(model::QuadraticMixing,x)
    return mixing_rule(model.func, x, model.prop,model.Aij)
end