const _vt = (
    VolumeAmount{MOLAR,VOLUME}(),
    Temperature(),
    SingleComponent())

const _vtx = (
    VolumeAmount{MOLAR,VOLUME}(),
    Temperature(),
    MaterialCompounds{MOLAR,FRACTION}())


function ∂a∂v(mt::SingleVT,model::HelmholtzModel,v,t)
    return ForwardDiff.derivative(dv->mol_helmholtz_impl(mt,model,dv,t),v)
    #return ForwardDiff.derivative(dv->mol_helmholtzR_impl(mt,model,dv,t),v) -RGAS*t/v
end

function ∂a∂v(mt::MultiVT,model::HelmholtzModel,v,t,x)
    return ForwardDiff.derivative(dv->mol_helmholtz_impl(mt,model,dv,t,x),v)
end

function ∂²a∂v²(mt::SingleVT,model::HelmholtzModel,v,t)
    return ForwardDiff.derivative(dv->∂a∂v(mt,model,dv,t),v)
    #return ForwardDiff.derivative(dv->mol_helmholtzR_impl(mt,model,dv,t),v) -RGAS*t/v
end

function ∂²a∂v²(mt::MultiVT,model::HelmholtzModel,v,t,x)
    return ForwardDiff.derivative(dv->∂a∂v(mt,model,dv,t,x),v)
end


function ∂a∂t(mt::SingleVT,model::HelmholtzModel,v,t)
    return ForwardDiff.derivative(dt->mol_helmholtz_impl(mt,model,v,dt),t)
end
    
function ∂a∂t(mt::MultiVT,model::HelmholtzModel,v,t,x)
    return ForwardDiff.derivative(dt->mol_helmholtz_impl(mt,model,v,dt,x),t)
end
    
function ∂a(mt::SingleVT,model::HelmholtzModel,v,t)
    a(x) = mol_helmholtz_impl(mt,model,first(x),last(x))
    T = promote_type(typeof(v),typeof(t))
    _v = T(v)
    _t = T(t)
    vt_vec =   SVector(_v,_t)
    ∂result = DiffResults.GradientResult(vt_vec)  
    _∂a =  ForwardDiff.gradient!(∂result, a,vt_vec)
    val_a =  DiffResults.value(_∂a)
    val_∂a = DiffResults.gradient(_∂a)
    return (val_∂a,val_a)
end

function ∂a(mt::MultiVT,model::HelmholtzModel,v,t,x)
    a(z) = mol_helmholtz_impl(mt,model,first(z),last(z),x)
    T = promote_type(typeof(v),typeof(t))
    _v = T(v)
    _t = T(t)
    vt_vec =   SVector(_v,_t)
    ∂result = DiffResults.GradientResult(vt_vec)  
    _∂a =  ForwardDiff.gradient!(∂result, a,vt_vec)
    val_a =  DiffResults.value(_∂a)
    val_∂a = DiffResults.gradient(_∂a)
    return (val_∂a,val_a)
end




function pressure_impl(mt::SingleVT,model::T,v,t) where T<: HelmholtzModel
    rho = one(v)/v
    dαdrho = ForwardDiff.derivative(drho->αR_impl(mt,model,drho,t),inv(v))
    return rho*(one(rho)+rho*dαdrho)*RGAS*t
end

function pressure_impl(mt::MultiVT,model::T,v,t,x) where T<: HelmholtzModel
    ρ = one(v)/v
    dαdρ = ForwardDiff.derivative(∂ρ->αR_impl(mt,model,∂ρ,t,x),inv(v))
    return ρ*(one(ρ)+ρ*dαdρ)*RGAS*t
end

function compressibility_factor_impl(mt::SingleVT,model::HelmholtzModel,v,t)
    p = pressure_impl(mt,model,v,t)
    return p*v/(RGAS*t)
end

function compressibility_factor_impl(mt::MultiVT,model::HelmholtzModel,v,t,x)
    p = pressure_impl(mt,model,v,t,x)
    return p*v/(RGAS*t)
end


function mol_entropy_impl(mt::SingleVT,model::HelmholtzModel,v,t)
    return -∂a∂t(mt,model,v,t)
end

function mol_entropy_impl(mt::MultiVT,model::HelmholtzModel,v,t,x)
    return -∂a∂t(mt,model,v,t,x)
end

function mol_enthalpy_impl(mt::SingleVT,model::HelmholtzModel,v,t)
    _da,aa =  ∂a(mt,model,v,t)
    _dv,_dt = _da
    return aa - _dv*v - _dt*t
end

function dpdt(mt::SingleVT,model::HelmholtzModel,v,t)
    p(z) = pressure_impl(mt,model,v,z)
    return ForwardDiff.derivative(p,t)
end

function mol_enthalpy_impl(mt::MultiVT,model::HelmholtzModel,v,t,x)
    _dv,_dt,aa =  ∂a(mt,model,v,t,x)
    return aa - _dv*v - dt*t
end

function mol_internal_energy_impl(mt::SingleVT,model::HelmholtzModel,v,t)
    _da,aa =  ∂a(mt,model,v,t)
    _dv,_dt = _da
    return aa  - dt*t
end

function mol_internal_energy_impl(mt::MultiVT,model::HelmholtzModel,v,t,x)
    _da,aa =  ∂a(mt,model,v,t,x)
    _dv,_dt = _da
    return aa  - dt*t
end

function mol_gibbs_impl(mt::SingleVT,model::HelmholtzModel,v,t)
    _da,aa =  ∂a(mt,model,v,t)
    _dv,_dt = _da
    return aa - _dv*v
end

function mol_gibbs_impl(mt::MultiVT,model::HelmholtzModel,v,t,x)
    _da,aa =  ∂a(mt,model,v,t,x)
    _dv,_dt = _da
    return aa - _dv*v
end

function ∂2a(mt::SingleVT,model::HelmholtzModel,v,t)
    a(x) = mol_helmholtz_impl(mt,model,first(x),last(x))
    T = promote_type(typeof(v),typeof(t))
    _v = T(v)
    _t = T(t)
    vt_vec =   SVector(_v,_t)
    ∂result = DiffResults.HessianResult(vt_vec)  
    _∂a =  ForwardDiff.hessian!(∂result, a,vt_vec)
    val_a =  DiffResults.value(_∂a)
    val_∂a = DiffResults.gradient(_∂a)
    val_∂2a = DiffResults.hessian(_∂a)
    return (val_∂2a,val_∂a,val_a)
end

function ∂2a(mt::MultiVT,model::HelmholtzModel,v,t,x)
    a(z) = mol_helmholtz_impl(mt,model,first(z),last(z),x)
    T = promote_type(typeof(v),typeof(t))
    _v = T(v)
    _t = T(t)
    vt_vec =   SVector(_v,_t)
    ∂result = DiffResults.HessianResult(vt_vec)
    _∂a =  ForwardDiff.hessian!(∂result, a,vt_vec)
    val_a =  DiffResults.value(_∂a)
    val_∂a = DiffResults.gradient(_∂a)
    val_∂2a = DiffResults.hessian(_∂a)
    return (val_∂2a,val_∂a,val_a)
end

function ∂2∂a2(mt::SingleVT,model::HelmholtzModel,v,t)
    a(x) = mol_helmholtz_impl(mt,model,first(x),last(x))
    T = promote_type(typeof(v),typeof(t))
    _v = T(v)
    _t = T(t)
    vt_vec =   SVector(_v,_t)
    return  ForwardDiff.hessian(a,vt_vec)
end


function ∂2∂a2(mt::MultiVT,model::HelmholtzModel,v,t,z)
    a(x) = mol_helmholtz_impl(mt,model,first(x),last(x),z)
    T = promote_type(typeof(v),typeof(t))
    _v = T(v)
    _t = T(t)
    vt_vec = SVector(_v,_t)
    return ForwardDiff.hessian(a,vt_vec)
end

function ∂p∂v(mt::SingleVT,model::HelmholtzModel,v,t)
return ForwardDiff.derivative(z->pressure_impl(mt,model,z,t),v)
end


function mol_isochoric_heat_capacity_impl(mt::SingleVT,model::HelmholtzModel,v,t)
    ∂2a = ∂2∂a2(mt,model,v,t)
    return -t*∂2a[2,2]
end

function mol_isochoric_heat_capacity_impl(mt::MultiVT,model::HelmholtzModel,v,t,x)
    ∂2a = ∂2∂a2(mt,model,v,t,x)
    return -t*∂2a[2,2]
end

function mol_isobaric_heat_capacity_impl(mt::SingleVT,model::HelmholtzModel,v,t)
    ∂2a = ∂2∂a2(mt,model,v,t)
    return t*(∂2a[1,2]^2/∂2a[1]-∂2a[2,2])
end

function mol_isobaric_heat_capacity_impl(mt::MultiVT,model::HelmholtzModel,v,t,x)
    ∂2a = ∂2∂a2(mt,model,v,t,x)
    return t*(∂2a[1,2]^2/∂2a[1]-∂2a[2,2])
end

function sound_speed_impl(mt::SingleVT,model::HelmholtzModel,v,t)
    ∂2a = ∂2∂a2(mt,model,v,t)
    Mr = molecular_weight(model)
    z = ThermoState.mw_div(one(MW),MW)
    return v*sqrt((∂2a[1]-∂2a[1,2]^2/∂2a[2,2])/Mr)
end


function sound_speed_impl(mt::MultiVT,model::HelmholtzModel,v,t,x)
    ∂2a = ∂2∂a2(mt,model,v,t,x)
    _MW = molecular_weight(model)
    Mr = mapreduce(ThermoState.mw_mul,+,x,_MW)
    return v*sqrt((∂2a[1]-∂2a[1,2]^2/∂2a[2,2])/Mr)
end


function isothermal_compressibility_impl(mt::SingleVT,model::HelmholtzModel,v,t)
    d2f = ∂²a∂v²(mt,model,v,t)
    return 1/v*d2f^-1
end

function isothermal_compressibility_impl(mt::MultiVT,model::HelmholtzModel,v,t,x)
    d2f = ∂²a∂v²(mt,model,v,t,x)
    return 1/v*d2f^-1
end

function pressure(model::HelmholtzModel,st::ThermodynamicState,unit = u"Pa")
    return pressure(state_type(st),model,st,unit)
end

function pressure(mt::SingleVT,model::HelmholtzModel,st::ThermodynamicState,unit = u"Pa")
    mw = molecular_weight(model)
    v = mol_volume(FromState(),st,u"m^3/mol",mw)
    t = temperature(FromState(),st)
    val = pressure_impl(mt,model,v,t)
    return convert_unit(u"Pa",unit,val)
end

function pressure(mt::MultiVT,model::HelmholtzModel,st::ThermodynamicState,unit = u"Pa")
    mw = molecular_weight(model)
    v = mol_volume(FromState(),st,u"m^3/mol",mw)
    t = temperature(FromState(),st)
    x = mol_fraction(FromState(),st,nothing,mw)
    val = pressure_impl(mt,model,v,t,x)
    return convert_unit(u"Pa",unit,val)
end

function mol_entropy(model::HelmholtzModel,st::ThermodynamicState,unit = u"J/(mol*K)")
    return mol_entropy(state_type(st),model,st,unit)
end

function mol_entropy(mt::SingleVT,model::HelmholtzModel,st::ThermodynamicState,unit = u"J/(mol*K)")
    mw = molecular_weight(model)
    v = mol_volume(FromState(),st,u"m^3/mol",mw)
    t = temperature(FromState(),st)
    val =  mol_entropy_impl(mt,model,v,t)
    return convert_unit(u"J/mol",unit,val)
end


function mol_entropy(mt::MultiVT,model::HelmholtzModel,st::ThermodynamicState,unit = u"J/(mol*K)")
    mw = molecular_weight(model)
    v = mol_volume(FromState(),st,u"m^3/mol",mw)
    t = temperature(FromState(),st,u"K")
    x = mol_fraction(FromState(),st,nothing,mw)
    val = mol_entropy_impl(mt,model,v,t,x)
    return convert_unit(u"J/(mol*K)",unit,val)
end

function mol_entropy(mt::SinglePT,model::HelmholtzModel,st::ThermodynamicState,unit = u"J/(mol*K)")
    mw = molecular_weight(model)
    p = pressure(FromState(),st)
    t = temperature(FromState(),st)
    phase = ThermoState.phase(FromState(),st)
    v = v_zero(mt,model,p,t,phase=phase)
    val =  mol_entropy_impl(mt,model,v,t)
    return convert_unit(u"J/mol",unit,val)
end

function mol_entropy(mt::MultiPT,model::HelmholtzModel,st::ThermodynamicState,unit = u"J/(mol*K)")
    mw = molecular_weight(model)
    p = pressure(FromState(),st)
    t = temperature(FromState(),st,u"K")
    phase = ThermoState.phase(FromState(),st)
    x = mol_fraction(FromState(),st,nothing,mw)
    v = v_zero(mt,model,p,t,x,phase=phase)
    val = mol_entropy_impl(mt,model,v,t,x)
    return convert_unit(u"J/(mol*K)",unit,val)
end

for mol_op in (:mol_helmholtz, :mol_gibbs, :mol_internal_energy, :mol_enthalpy)
    mol_op_impl = Symbol(mol_op,:_impl)
    @eval begin 
        function $mol_op_impl(model::HelmholtzModel,args...)
            return $mol_op_impl(model_type(model),model,args...)
        end

        function $mol_op(model::HelmholtzModel,st::ThermodynamicState,unit=u"J/mol")
            return $mol_op(state_type(st),model,st,unit)
        end    

        function $mol_op(mt::SingleVT,model::HelmholtzModel,st::ThermodynamicState,unit = u"J/mol")
            mw = molecular_weight(model)
            v = mol_volume(FromState(),st,u"m^3/mol",mw)
            t = temperature(FromState(),st)
            val = $mol_op_impl(mt,model,v,t)
            return convert_unit(u"J/mol",unit,val)
        end

        function $mol_op(mt::MultiVT,model::HelmholtzModel,st::ThermodynamicState,unit = u"J/mol")
            mw = molecular_weight(model)
            v = mol_volume(FromState(),st,u"m^3/mol",mw)
            t = temperature(FromState(),st)
            x = mol_fraction(FromState(),st,nothing,mw)
            val = $mol_op_impl(mt,model,v,t,x)
            return convert_unit(u"J/mol",unit,val)
        end

        function $mol_op(mt::SinglePT,model::HelmholtzModel,st::ThermodynamicState,unit = u"J/mol")
            mw = molecular_weight(model)
            p = pressure(FromState(),st)
            t = temperature(FromState(),st)
            phase = ThermoState.phase(FromState(),st)
            v = v_zero(mt,model,p,t,phase=phase)           
            val = $mol_op_impl(mt,model,v,t)
            return convert_unit(u"J/mol",unit,val)
        end

        function $mol_op(mt::MultiPT,model::HelmholtzModel,st::ThermodynamicState,unit = u"J/mol")
            mw = molecular_weight(model)
            p = pressure(FromState(),st)
            t = temperature(FromState(),st)
            x = mol_fraction(FromState(),st,nothing,mw)
            phase = ThermoState.phase(FromState(),st)
            v = v_zero(mt,model,p,t,x,phase=phase)
            val = $mol_op_impl(mt,model,v,t,x)
            return convert_unit(u"J/mol",unit,val)
        end
    end
end

for (mol_op,mass_op) in zip((:mol_helmholtz, :mol_gibbs, :mol_internal_energy, :mol_enthalpy),(:mass_helmholtz, :mass_gibbs, :mass_internal_energy, :mass_enthalpy))
    @eval begin 
        function $mass_op(model::HelmholtzModel,st::ThermodynamicState,unit=u"J/kg")
            mw = molecular_weight(model)
            mol_val = $mol_op(model,st)
            val = mol_val/molar_mass(FromState(),st,u"kg/mol",mw)
            return convert_unit(u"J/kg",unit,val)
        end    
    end
end

for (mol_op,total_op) in zip((:mol_helmholtz, :mol_gibbs, :mol_internal_energy, :mol_enthalpy),(:total_helmholtz, :total_gibbs, :total_internal_energy, :total_enthalpy))
    @eval begin 
        function $total_op(model::HelmholtzModel,st::ThermodynamicState,unit=u"J")
            mw = molecular_weight(model)
            mol_val = $mol_op(model,st)
            val = mol_val*moles(FromState(),st,u"mol",mw)
            return convert_unit(u"J",unit,val)
        end    
    end
end


function total_entropy(model::HelmholtzModel,st::ThermodynamicState,unit=u"J/(kg*K)")
    mw = molecular_weight(model)
    mol_val = mol_entropy(model,st)
    val = mol_val*moles(FromState(),st,u"mol",mw)
    return convert_unit(u"J/(kg*K)",unit,val)
end 

function mol_volume(model::HelmholtzModel,st::ThermodynamicState,unit = u"m^3/mol")
    return mol_volume(state_type(st),model,st,unit)
end

function mol_volume(mt::SinglePT,model::HelmholtzModel,st::ThermodynamicState,unit = u"m^3/mol")
    mw = molecular_weight(model)
    p = pressure(FromState(),st)
    t = temperature(FromState(),st)
    phase = ThermoState.phase(FromState(),st)
    v = v_zero(volume_solver_type(model),mt,model,p,t,phase=phase)           
    return convert_unit(u"m^3/mol",unit,v)
end


function mol_volume(mt::MultiPT,model::HelmholtzModel,st::ThermodynamicState,unit = u"m^3/mol")
    mw = molecular_weight(model)
    p = pressure(FromState(),st)
    t = temperature(FromState(),st)
    x = mol_fraction(FromState(),st,nothing,mw)
    phase = ThermoState.phase(FromState(),st)
    v = v_zero(volume_solver_type(model),mt,model,p,t,x,phase=phase)           
    return convert_unit(u"m^3/mol",unit,v)
end


function mol_density(model::HelmholtzModel,st::ThermodynamicState,unit =u"mol/(m^3)")
    v = mol_volume(model,st)
    res =  one(v)/v
    return convert_unit(u"mol/(m^3)",unit,res)
end

function total_volume(model::HelmholtzModel,st::ThermodynamicState,unit = u"m^3")
    mw = molecular_weight(model)
    v  = mol_volume(model,st)
    n = moles(FromState(),st,u"mol",mw)
    res =  v*n
    return convert_unit(u"m^3",unit,res)
end

function mass_volume(model::HelmholtzModel,st::ThermodynamicState,unit = u"m^3/kg")
    mw = molecular_weight(model)
    v  = mol_volume(model,st)
    val = v/molar_mass(FromState(),st,u"kg/mol",mw)
    return convert_unit(u"m^3/kg",unit,val)
end

function mass_density(model::HelmholtzModel,st::ThermodynamicState,unit =u"kg/(m^3)")
    v  = mass_volume(model,st)
    res =  one(v)/v
    return convert_unit(u"kg/(m^3)",unit,res)
end

function compressibility_factor(model::HelmholtzModel,st::ThermodynamicState)
    return compressibility_factor(state_type(st),model,st)
end

function compressibility_factor(mt::SingleVT,model::HelmholtzModel,st::ThermodynamicState)
    v =mol_volume(FromState(),st)
    t = temperature(FromState(),st)
    p = pressure_impl(mt,model,v,t)
    return p*v/(RGAS*t)
end

function compressibility_factor(mt::MultiVT,model::HelmholtzModel,st::ThermodynamicState)
    v =mol_volume(FromState(),st)
    t = temperature(FromState(),st)
    x = mol_fraction(FromState(),st)
    p = pressure_impl(mt,model,v,t,x)
    return p*v/(RGAS*t)
end