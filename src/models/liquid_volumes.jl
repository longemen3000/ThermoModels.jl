struct Rackett{MW,T} <: SatLiquidVolumeModel
    mw::MW
    tc::T
    pc::T
    zc::T
    
end


Rackett(;tc,pc,zc,mw=nothing) = Rackett(tc,pc,zc,mw)


function Rackett(model::ThermoModel)
    mw = molecular_weight(model)
    tc = temperature(model,CriticalPoint())
    pc = pressure(model,CriticalPoint())
    zc = compressibility_factor(model,CriticalPoint())
    return Rackett(tc,pc,zc,mw)
end

function mol_volume_impl(::SingleSatT,model::Rackett,t::T) where T
    _2_7 = T(2//7)
    _1 = T(1.0)
    tc = only(model.tc)
    zc = only(model.zc)
    pc = only(model.pc)
    RGAS*tc/pc*zc^(_1 + (_1 - T/Tc)^(_2_7))
end


function mol_volume_impl(::MultiSatT,model::Rackett,t::T,x) where T
    tc = model.tc
    zc = model.zc
    pc = model.pc
    _2_7 = T(2//7)
    _1 = T(1.0)
    zmix = dot(zc,x)
    tmix = dot(tc,x)
    mw = dot(MW,x)
    tr = t/tmix
    RGAS*tc/pc*zmix^(_1 + (_1 - tr)^(_2_7))
end





