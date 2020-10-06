#volume_solver(model::HelmholtzModel,args...;kwargs...) = volume_solver(model_type(model),volume_solver_type(model),model,args...;kwargs...)

function volume_solver(
    mt::SinglePT,
    method::VolumeBisection,
    model::HelmholtzModel,
    p,
    t,
    v0 = nothing;
    no_pts = 7)

    _p(z) = pressure_impl(QuickStates.vt(),model, z, t)
    fp(z) = _p(z) - p
    dfp(z) = ForwardDiff.derivative(fp, z)
    if isnothing(v0)
        min_v = only(covolumes(model))
        max_v = 40 * t / p #approx 5 times ideal gas
    #this is to be sure that (min_v,max_v) is a bracketing interval
        while fp(max_v) > 0
            max_v *= 2
            (max_v > typemax(max_v)/8) && break
        end
        while fp(min_v) < 0
            min_v *= 0.5
            (min_v < sqrt(eps((max_v)))) && break
        end
        vv = find_zeros(fp, min_v, max_v, no_pts = no_pts)
        return (first(vv),last(vv))
    else
        try
            vv =  Roots.find_zero(fp, v0)
            return (vv,vv)
        catch
            return volume_solver(mt,method,model, p, t; no_pts = no_pts)
        end
    end
end

function volume_solver(
    mt::MultiPT,
    method::VolumeBisection,
    model::HelmholtzModel,
    p,
    t,
    x,
    v0 = nothing;
    no_pts = 7)
    _p(z) = pressure_impl(QuickStates.vtx(),model, z, t, x)
    fp(z) = _p(z) - p
    dfp(z) = ForwardDiff.derivative(fp, z)
    if isnothing(v0)
        min_v = dot(covolumes(model),x)
        max_v = 40 * t / p #approx 5 times ideal gas
    #this is to be sure that (min_v,max_v) is a bracketing interval
        while fp(max_v) > 0
            max_v *= 2
            (max_v > typemax(max_v)/8) && break
        end
        while fp(min_v) < 0
            min_v *= 0.5
            (min_v < sqrt(eps((max_v)))) && break
        end
        vv = find_zeros(fp, min_v, max_v, no_pts = no_pts)
        return (first(vv),last(vv))
    else
        try
            vv =  find_zero(fp, v0)
            return (vv,vv)
        catch
            return volume_solver(mt,method,model, p, t,x; no_pts = no_pts)
        end
    end
end


function flash_impl(mt::SingleSatP,model::HelmholtzModel, _p)

    p = float(normalize_units(_p))
    TTT = typeof(p)
    function v_solver(p0,t0,v0=nothing,pts=7) 
        return volume_solver(QuickStates.pt(),volume_solver_type(model),model,p0,t0,v0,no_pts=pts)
    end
    if (Pc = pressure(model,CriticalPoint())) ≈ p
        vc = mol_volume(model,CriticalPoint())
        tc=temperature(model,CriticalPoint())
        return state(mol_v=vc,t=tc)
    elseif p > Pc
        throw(error("the phase is supercritical at pressure = $p Pa"))
    end
    Tc = temperature(model,CriticalPoint())
    model_pred = SingleSatPredictor(model)
    t_pred = temperature_impl(mt,model_pred,p)
    vv = v_solver(p,t_pred,nothing,21)
    if length(vv) <= 1
        throw(error("the phase is stable at pressure = $p Pa"))
    elseif length(vv) >= 2
        v1 = first(vv)
        v2 = last(vv)
        px = p
        _A(z, t) = mol_helmholtz_impl(QuickStates.vt(),model, z, t)
        v1old = 0.0
        v2old = 0.0
        Told = -100.0
        Tx = t_pred
        for i = 1:20
            if i > 1
                v1old = v1
                v2old = v2
            end
            A = z -> _A(z, Tx)
            _v1 = Threads.@spawn v_solver($p, $Tx, 0.9 * $v1)::Tuple{TTT,TTT}
            _v2 = Threads.@spawn v_solver($p, $Tx, 1.1 * $v2)::Tuple{TTT,TTT}
            v1 = first(fetch(_v1))
            v2 = last(fetch(_v2))
            if abs(v1 - v1old) / v1 < 1e-15 && i > 1
                #println("v1 condition")
                break
            elseif abs(v2 - v2old) / v2 < 1e-15 && i > 1
                #println("v2 condition")
                break
            end
            
            _px(T) = (_A(v1, T) - _A(v2, T)) / (v2 - v1) - p
            #_px(T) = exp(_A(v1,T)-_A(v2,T)) - one(Tx)
            Told = Tx
            Tx = Roots.find_zero(_px, Tx)
            if abs(Told - Tx) < 1e-10 * Tx
                #println("P condition")
                break
            end
        end

        phase1 = state(mol_v = v1, t=Tx ,phase=:liquid)
        phase2 = state(mol_v = v2, t=Tx, phase=:gas)
        return (phase1, phase2)
    end
end

function flash_impl(mt::SingleSatT,model::HelmholtzModel, _t)
    t = float(normalize_units(_t))
    tc = ThermoState.temperature(model,CriticalPoint())
    function v_solver(p0,t0,v0=nothing,pts=7) 
        return volume_solver(QuickStates.pt(),volume_solver_type(model),model,p0,t0,v0,no_pts=pts)
    end
    if tc ≈ t
        return state(mol_v=mol_volume(model,CriticalPoint()),t=tc)
    elseif t > tc
        throw(error("the phase is supercritical at temperature = $T K"))
    end
    model_pred = SingleSatPredictor(model)
    p_pred = pressure_impl(mt,model_pred,t)
    p = p_pred
    _p(z) = pressure_impl(QuickStates.vt(),model, z, t)
    fp(z) =_p(z) - p
    min_v = only(covolumes(model))
    max_v = 40 * t / p_pred #approx 5 times ideal gas
    #this is to be sure that (min_v,max_v) is a bracketing interval
    while fp(max_v) > 0
        max_v *= 2
    end
    while fp(min_v) < 0
        min_v *= 0.5
    end
    vv = Roots.find_zeros(fp, min_v, max_v, no_pts = 21)
    if vv[begin] ≈ vv[end]
        throw(error("the phase is stable at temperature = $T K"))
    else
        v1 = vv[begin]
        v2 = vv[end]
        p1 = _p(v1)
        p2 = _p(v2)
        px = p_pred
        A(z) = fugacity_coeff_impl(QuickStates.vt(),model, z, t)
        v1old = 0
        v2old = 0
        for i = 1:20
            if i > 1
                v1old = v1
                v2old = v2
            end
            _v1 = Threads.@spawn v_solver($px, $t, 0.9 * $v1)
            _v2 = Threads.@spawn v_solver($px, $t, 1.1 * $v2)

            v1 = first(fetch(_v1))
            v2 = last(fetch(_v2))
            if abs(v1 - v1old) / v1 < 1e-15 && i > 1
                #println("v1 condition")
                break
            elseif abs(v2 - v2old) / v2 < 1e-15 && i > 1
                #println("v2 condition")
                break
            end
            pold = px
            px = px*exp(A(v1)-A(v2))
            #px = (A(v1) - A(v2)) / (v2 - v1)
            #@show px
            if abs(pold - px) < 1e-10 * px
                #println("P condition")
                break
            end

        end
        #exp(logphi)*P= f
        phase1 = state(mol_v = v1, t= t,phase=:liquid)
        phase2 = state(mol_v = v2, t= t,phase=:gas)
        return (phase1, phase2)
    end
end

function equilibria(model::HelmholtzModel,st::ThermodynamicState)
    return equilibria(state_type(st),model,st)
end

function equilibria(mt::SingleSatP,model::HelmholtzModel,st::ThermodynamicState)
    p = pressure(FromState(),st)
    return flash_impl(mt,model,p)
end

function equilibria(mt::SingleΦP,model::HelmholtzModel,st::ThermodynamicState)
    p = pressure(FromState(),st)
    ϕ = quality(FromState(),st)
    m = mass(FromState()st,u"kg",molecular_weight(model))
    mv = ϕ*m
    ml = (1-ϕ)*m
    st_l, st_g = flash_impl(mt,model,p)
    Teq = temperature(FromState(),st_l)
    vl = mol_volume(FromState(),st_l)
    vg = mol_volume(FromState(),st_g)
    tagl = :liquid
    tagv = :gas
    st1 = state(mass=ml,t=Teq,mol_v=vl,phase=tagl)
    st2 = state(mass=mv,t=Teq,mol_v=vv,phase=tagv)
    return (st1,st2)
end

function equilibria(mt::SingleSatT,model::HelmholtzModel,st::ThermodynamicState)
    t = temperature(FromState(),st)
    return flash_impl(mt,model,t)
end

function equilibria(mt::SingleΦT,model::HelmholtzModel,st::ThermodynamicState)
    t = temperature(FromState(),st)
    ϕ = quality(FromState(),st)
    m = mass(FromState()st,u"kg",molecular_weight(model))
    mv = ϕ*m
    ml = (1-ϕ)*m
    st_l, st_g = flash_impl(mt,model,t)
    Peq = pressure(FromState(),st_l)
    vl = mol_volume(FromState(),st_l)
    vg = mol_volume(FromState(),st_g)
    tagl = :liquid
    tagv = :gas
    st1 = state(mass=ml,p=Peq,mol_v=vl,phase=tagl)
    st2 = state(mass=mv,p=Peq,mol_v=vv,phase=tagv)
    return (st1,st2)
end

function pressure(mt::SingleSatT,model::HelmholtzModel,st::ThermodynamicState,unit)
    t = temperature(FromState(),st)
    st_l, st_g = flash_impl(mt,model,t)
    p1 = pressure(model,st_l)
    p2 = pressure(model,st_g)
    Peq = 0.5*(p1+p2)
    return convert_unit(u"Pa",unit,Peq)
end

function temperature(mt::SingleSatP,model::HelmholtzModel,st::ThermodynamicState,unit)
    p = pressure(FromState(),st)
    st_l, st_g = flash_impl(mt,model,p)
    Teq = temperature(FromState(),st_l)
    return convert_unit(u"K",unit,Teq)
end

function pressure(mt::SingleΦT,model::HelmholtzModel,st::ThermodynamicState,unit)
    return pressure(QuickStates.sat_t(),model,st,unit)
end

function temperature(mt::SingleΦP,model::HelmholtzModel,st::ThermodynamicState,unit)
    return temperature(QuickStates.sat_p(),model,st,unit)
end