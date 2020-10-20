function volume_solver(
    mt::SinglePT,
    method::VolumeBisection,
    model::HelmholtzModel,
    p,
    t,
    v0 = nothing;
    phase=:liquid)
    no_pts = method.pts
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
            if phase == :liquid
            vv =  Roots.find_zero(fp, v0)
            else
                vv = gas_fzero(_p,p,v0,t)
                if isnan(vv)
                    println(v0)
                    return volume_solver(mt,method,model, p, t; no_pts = no_pts)
                end
            end
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



#returns a mol volume from P and T
function v_zero(mt::SinglePT,model::HelmholtzModel,p,t,v0=nothing;phase=:unspecified)
    _vt = QuickStates.vt()
    _p(_v) = pressure_impl(_vt,model,_v,t)
    _g(_v) = mol_gibbs_impl(_vt,model,_v,t)
    if is_gas(phase)
        if isnothing(v0)
            v = gas_fzero(_p,p,t)
        else
            v = gas_fzero(_p,p,t,v0) 
        end
        if isnan(v)
            #println("isnan with v=",v0,", p= ",p,", t = ",t)
            vs = v_zeros(mt,volume_solver_type(model),model,p,t)
            return last(vs)
        else
            return v
        end
    elseif is_liquid(phase)
        if isnothing(v0)
            vs = v_zeros(mt,volume_solver_type(model),model,p,t)
            v = first(vs)
            return v
        else
            try
                v= roots_fzero(_p,p,t,v0)
                return v
            catch
                println("isnan with v=",v0,", p= ",p,", t = ",t)
                vs = v_zeros(mt,volume_solver_type(model),model,p,t)
                v = first(vs)
                return v
            end
        end
    else #return the phase with least gibbs energy
        if isnothing(v0)
            vs = v_zeros(mt,volume_solver_type(model),model,p,t)
            v1 = first(vs)
            v2 = last(vs)
            (v1 == v2) && return v1
            g1 = _g(v1)
            g2 = _g(v2)
            if g1 < g2
                return v1
            elseif g2 < g1
                return v2
            else #different volumes, but equal energies, equilibria conditions
                return v2
            end
        else
            try
                return roots_fzero(_p,p,t,v0)
            catch
                vs = v_zeros(mt,volume_solver_type(model),model,p,t)
                v1 = first(vs)
                v2 = last(vs)
                (v1 == v2) && return v1
                g1 = _g(v1)
                g2 = _g(v2)
                if g1 < g2
                    return v1
                elseif g2 < g1
                    return v2
                else #different volumes, but equal energies, equilibria conditions
                    return v2
                end
            end
        end
    end
end
#returns all mol volumes from P and T
function v_zeros(::SinglePT,
    method::VolumeBisection,
    model::HelmholtzModel,
    p,
    t)
    no_pts = method.pts
    _p(z) = pressure_impl(QuickStates.vt(),model, z, t)
    fp(z) = _p(z) - p
    dfp(z) = ForwardDiff.derivative(fp, z)
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
end

function v_zero(mt::MultiPT,model::HelmholtzModel,p,t,x,v0=nothing;phase=:unspecified)
    _vt = QuickStates.vtx()
    _p(_v) = pressure_impl(_vt,model,_v,t,x)
    _g(_v) = mol_gibbs_impl(_vt,model,_v,t,x)
    if is_gas(phase)
        if isnothing(v0)
            v = gas_fzero(_p,p,t)
        else
            v = gas_fzero(_p,p,t,v0) 
        end
        
        if isnan(v)
            #println("isnan with v=",v0,", p= ",p,", t = ",t)
            vs = v_zeros(mt,volume_solver_type(model),model,p,t,x)
            return last(vs)
        else
            return v
        end

    elseif is_liquid(phase)
        if isnothing(v0)
            vs = v_zeros(mt,volume_solver_type(model),model,p,t,x)
            v = first(vs)
            return v
        else
            try
                v= roots_fzero(_p,p,t,v0)
                return v
            catch
                vs = v_zeros(mt,volume_solver_type(model),model,p,t,x)
                v = first(vs)
                return v
            end
        end
    else #return the phase with least gibbs energy
        if isnothing(v0)
            vs = v_zeros(mt,volume_solver_type(model),model,p,t,x)
            v1 = first(vs)
            v2 = last(vs)
            (v1 == v2) && return v1
            g1 = _g(v1)
            g2 = _g(v2)
            if g1 < g2
                return v1
            elseif g2 < g1
                return v2
            else #different volumes, but equal energies, equilibria conditions
                return v2
            end
        else
            try
                return roots_fzero(_p,p,t,v0)
            catch
                vs = v_zeros(mt,volume_solver_type(model),model,p,t,x)
                v1 = first(vs)
                v2 = last(vs)
                (v1 == v2) && return v1
                g1 = _g(v1)
                g2 = _g(v2)
                if g1 < g2
                    return v1
                elseif g2 < g1
                    return v2
                else #different volumes, but equal energies, equilibria conditions
                    return v2
                end
            end
        end
    end
end
#returns all mol volumes from P and T
function v_zeros(::MultiPT,
    method::VolumeBisection,
    model::HelmholtzModel,
    p,
    t,
    x)
    no_pts = method.pts
    _p(z) = pressure_impl(QuickStates.vtx(),model, z, t,x)
    fp(z) = _p(z) - p
    dfp(z) = ForwardDiff.derivative(fp, z)
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
    v1 = first(vv)
    v2 = last(vv)
    #v1 =  Roots.find_zero(fp,first(vv))
    #v2 = Roots.find_zero(fp,last(vv))
    return (v1,v2)
end

function roots_fzero(p,pspec,t,v)
    f0(x) = p(x) - pspec
    return Roots.find_zero(f0,v)
end

function gas_fzero(p,Pspec,t,v0 = 40 * t / Pspec;
    relax=0.9*one(float(Pspec)),
    maxevals=100,
    atol = zero(float(Pspec)),
    rtol = 8*eps(one(float(Pspec))))

    _Pspec,_v0,_t,_atol,_rtol,_relax = promote(Pspec,v0,t,atol,rtol,relax)
    return _gas_fzero(p,_Pspec,_v0,_t,_relax,maxevals,_atol,_rtol)
end

function _gas_fzero(P,Pspec::TT,v0::TT,T::TT,relax = TT(0.9),maxevals=100, atol=zero(TT), rtol=8*eps(TT)) where TT
    pspec = Pspec
    R = RGAS #constant defined outside, 8.314...
    v₀ = v0
    P₀= P(v0)
    
    λ = relax
    _1 = one(TT)
    _0 = zero(TT)
    P₁ = _0
    v₁ = _0
    


    nan = _0/_0
    Pold = nan
    vold = nan

    RT = R*T
    RT_P = RT/Pspec
    invRT = _1/RT
    preA = Pspec*invRT*invRT
    preB = Pspec*invRT
    count = 0
    uatol = atol / oneunit(atol) * oneunit(real(v₁))
    adjustunit = oneunit(real(P₀))/oneunit(real(v₀))


    while count < maxevals
        #at high volumes, P ≈ RT/v - a/v2
        a = v₀*(RT-P₀*v₀)
        Δ = RT_P^2 - 4*a/Pspec
        if Δ > 0
            v₁ = 0.5*(RT_P + sqrt(Δ))
            P₁ = P(v₁)
        else
            if !isnan(vold)
            v₁ = 0.5*(v₀+vold)
            P₁ = P(v₁)
            else
                v₁ = 20v₀
                P₁ = P(v₁)
            end
        end

        if (P₁ > Pspec)
            vmax = v₁
            Pmax = P₁
        end
        if (Pspec < P₀ < P₁) #step in the wrong direction, reduce volume directly
            v₁ = v₀*1.1
            P₁ = P(v₁)
        end
        ΔP₀ = Pspec - P₀
        ΔP = Pspec - P₁
        ΔPold = Pspec - Pold

        #@show ΔP , ΔP₀, ΔPold,v₁
        abs(ΔP) <= eps(Pspec) && return v₁
        abs(ΔP) <= adjustunit * max(uatol, abs(v₁) * rtol) && return v₁
        iszero(ΔP) && return v₁
        isnan(ΔP) && return nan
        if (ΔP == ΔP₀)
            println("weird, v₁ =",v₁," and v₀ =",v₀)
            sdp = sign(ΔP)
            v_sprev = sign(P(prevfloat(v₁))-Pspec)
            v_snext = sign(P(prevfloat(v₁))-Pspec)
            p_sprev = sign(prevfloat(ΔP))
            p_snext = sign(nextfloat(ΔP))
            sdp * v_sprev <= 0 && return v₁
            sdp * v_snext <= 0 && return v₁
            sdp * p_sprev <= 0 && return v₁
            sdp * p_snext <= 0 && return v₁
            return nan
        end
        Pold,vold,P₀,v₀ = P₀,v₀,P₁,v₁
        count +=1 
    end
    return nan
end
