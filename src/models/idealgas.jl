const SingleAny = Tuple{T1,T2,SingleComponent} where {T1,T2}
const MultiAny = Tuple{T1,T2,MaterialCompounds} where {T1,T2}

struct IdealGas{T,MW,CP} <: GasModel
    type::T 
    mw::MW 
    cp::CP
end

IdealGas() = IdealGas(SingleComponent(),nothing,nothing)
function IdealGas(mw)
    if length(mw) == 1
        return IdealGas(SingleComponent(),only(mw),nothing)
    else
        return IdealGas(MaterialCompounds{MOLAR,FRACTION}(),mw,nothing)
    end
end

function IdealGas(mw,cp)
    if length(mw) == 1
        return IdealGas(SingleComponent(),only(mw),cp)
    else
        return IdealGas(MaterialCompounds{MOLAR,FRACTION}(),mw,cp)
    end
end

function pressure_impl(mt::SingleVT,model::IdealGas,v,t)
    return RGAS*t/v
end

function pressure_impl(mt::MultiVT,model::IdealGas,v,t,x)
    return RGAS*t/v
end

function temperature_impl(mt::SinglePV,model::IdealGas,p,v)
    return RGAS*t/v
end


function temperature_impl(mt::MultiPV,model::IdealGas,p,v,x)
    return p*v/RGAS
end

function compressibility_factor(mt::SingleAny,model::IdealGas,a1,a2)
    return one(a1)*one(a2)
end

function compressibility_factor(mt::MultiAny,model::IdealGas,a1,a2,a3)
    return one(a1)*one(a2)*one(eltype(a3))
end


function mol_volume_impl(mt::SinglePT,model::IdealGas,p,t)
    return RGAS*T/p
end

function mol_volume_impl(mt::MultiPT,model::IdealGas,p,t,x)
    return RGAS*T/p
end









