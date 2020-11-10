#returns a mol volume from P and T
function v_zero(mt::SinglePT,model,args...;kwargs...)
    return v_zero(volume_solver_type(model),mt,model,args...;kwargs...)
end

function v_zero(mt::MultiPT,model,args...;kwargs...)
    return v_zero(volume_solver_type(model),mt,model,args...;kwargs...)
end

function v_zero(::SUVA,mt::SinglePT,model,p,t,v0=nothing;phase=:unspecified)
    _vt = QuickStates.vt()
    _p(_v) = pressure_impl(_vt,model,_v,t)
    if is_gas(phase)
        if isnothing(v0)
            v = gas_fzero(_p,p,t)
        else
            v = gas_fzero(_p,p,t,v0) 
        end
    elseif is_liquid(phase)
        if isnothing(v0)
            v = v_zero(VolumeBisection(),mt,model,p,t;phase=phase)
            #v = liquid_fzero(_p,p,t)
            return v
        else
            try
                v= roots_fzero(_p,p,t,v0)
                return v
            catch
                v = v_zero(VolumeBisection(),mt,model,p,t;phase=phase)
                #v = liquid_fzero(_p,p,t)
                return v
            end
        end

    else #return the phase with least gibbs energy
        if isnothing(v0)
            return v_zero_general(mt,model,p,t)
        else
            try
                return roots_fzero(_p,p,t,v0)
            catch
                return v_zero_general(mt,model,p,t)
            end
        end
    end
end
#returns all mol volumes from P and T

function v_zero(::SUVA,mt::MultiPT,model,p,t,x,v0=nothing;phase=:unspecified)
    _vt = QuickStates.vtx()
    _p(_v) = pressure_impl(_vt,model,_v,t,x)
    if is_gas(phase)
        if isnothing(v0)
            v = gas_fzero(_p,p,t)
        else
            v = gas_fzero(_p,p,t,v0) 
        end

    elseif is_liquid(phase)
        if isnothing(v0)
            v = v_zero(VolumeBisection(),mt,model,p,t,x;phase=phase)
            #v = liquid_fzero(_p,p,t)
            return v
        else
            try
                v= roots_fzero(_p,p,t,v0)
                return v
            catch
                v = v_zero(VolumeBisection(),mt,model,p,t,x;phase=phase)
                #v = liquid_fzero(_p,p,t)
                return v
            end
        end
    else #return the phase with least gibbs energy
        if isnothing(v0)
            return v_zero_general(mt,model,p,t,x)

        else
            try
                return roots_fzero(_p,p,t,v0)
            catch

                return v_zero_general(mt,model,p,t,x)
            end
        end
    end
end
#returns all mol volumes from P and T
function v_zeros(mt,model,args...;kwargs...)
    return v_zeros(volume_solver_type(model),mt,model,args...;kwargs...)
end
function v_zeros(::SUVA
    ,mt::SinglePT
    ,model::HelmholtzModel
    ,p
    ,t)

    _vt = QuickStates.vt()
    _p(_v) = pressure_impl(_vt,model,_v,t)
    v1 = liquid_fzero(_p,p,t)
    v2 = gas_fzero(_p,p,t)
    return (v1,v2)
end

function v_zeros(::SUVA
    ,mt::MultiPT
    ,model
    ,p
    ,t
    ,x)

    _vt = QuickStates.vtx()
    _p(_v) = pressure_impl(_vt,model,_v,t,x)
    v1 = liquid_fzero(_p,p,t)
    v2 = gas_fzero(_p,p,t)
    return (v1,v2)
end


function v_zeros(vmethod::VolumeBisection
    ,mt::SinglePT
    ,model
    ,p
    ,t
    ;phase=:unspecified)
    no_pts = vmethod.pts
    _p(z) = pressure_impl(QuickStates.vt(),model, z, t)
    fp(z) = _p(z) - p
        max_v = 40 * t / p #approx 5 times ideal gas
    #this is to be sure that (min_v,max_v) is a bracketing interval
    while (fp(max_v)) > 0 |  (_p(max_v) < 0)
        max_v *= 2
        (max_v > typemax(max_v)*0.125) && break
    end

    vv = Roots.find_zeros(fp, eps(typeof(max_v)), max_v, no_pts = no_pts)
    return (first(vv),last(vv))
end



function v_zeros(vmethod::VolumeBisection,
    mt::MultiPT,
    model,
    p,
    t,
    x
    ;phase=:unspecified)
    no_pts = vmethod.pts
    _p(z) = pressure_impl(QuickStates.vtx(),model, z, t,x)
    fp(z) = _p(z) - p
        max_v = 40 * t / p #approx 5 times ideal gas
    #this is to be sure that (min_v,max_v) is a bracketing interval
    while (fp(max_v)) > 0 |  (_p(max_v) < 0)
        max_v *= 2
        (max_v > typemax(max_v)*0.125) && break
    end
 
    vv = Roots.find_zeros(fp, eps(typeof(max_v)), max_v, no_pts = no_pts)
    v1 = first(vv)
    v2 = last(vv)
    #v1 =  Roots.find_zero(fp,first(vv))
    #v2 = Roots.find_zero(fp,last(vv))
    return (v1,v2)
end

function v_zero(vmethod::VolumeBisection
    ,mt::SinglePT
    ,model
    ,p
    ,t
    ,v0=nothing
    ;phase=:unspecified)
    
    _p(z) = pressure_impl(QuickStates.vt(),model, z, t)

    if !(v0 === nothing)
        try
            return roots_fzero(_p,p,t,v0)
        catch
            return v_zero_general(mt,model,p,t)
        end
    end
    gas = is_gas(phase)
    liquid = is_liquid(phase)
    if !gas & !liquid
        return v_zero_general(mt,model,p,t)
    end
       
    (vl,vg) = v_zeros(vmethod,mt,model,p,t)
    if gas
        return vg
    else
        return vl
    end
end


function v_zero(vmethod::VolumeBisection
    ,mt::MultiPT
    ,model
    ,p
    ,t
    ,x
    ,v0=nothing
    ;phase=:unspecified)
    
    _p(z) = pressure_impl(QuickStates.vtx(),model, z, t,x)

    if !(v0 === nothing)
        try
            return roots_fzero(_p,p,t,v0)
        catch
            return v_zero_general(mt,model,p,t,x)
        end
    end
    gas = is_gas(phase)
    liquid = is_liquid(phase)
    if !gas & !liquid
        return v_zero_general(mt,model,p,t,x)
    end
       
    (vl,vg) = v_zeros(vmethod,mt,model,p,t,x)
    if gas
        return vg
    else
        return vl
    end
end



function v_zero_general(mt::SinglePT,model,p,t)
    _vt = QuickStates.vt()
    _g(_v) = mol_gibbs_impl(_vt,model,_v,t)
    _p(_v) = pressure_impl(_vt,model,_v,t)  
    
    v1,v2 = v_zeros(mt,model,p,t)
    #(v1 == v2) && return v1
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

#finds the root with lowest gibbs energy, corrected for gas
function v_zero_general(mt::MultiPT,model,p,t,x)
    _vt = QuickStates.vtx()
    _g(_v) = mol_gibbs_impl(_vt,model,_v,t,x)
    _p(_v) = pressure_impl(_vt,model,_v,t,x)  
    
    v1,v2 = v_zeros(mt,model,p,t,x)
    #(v1 == v2) && return v1
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

function v_zero_general(mt::MultiPT,model::CubicModel,p,t,x)
    _vt = QuickStates.vtx()
    _g(_v) = αR_impl(_vt,model,_v,t,x) - p*_v/(RGAS*t)
    _p(_v) = pressure_impl(_vt,model,_v,t,x)  
    
    v1,v2 = v_zeros(mt,model,p,t,x)
    #(v1 == v2) && return v1
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

function roots_fzero(p,pspec,t,v)
    f0(x) = p(x) - pspec
    return Roots.find_zero(f0,v)
end

#Sucessive van der waals approximations, patent pending (lol)
#aproximates locally the function to a van der waals
#gas_fzero aproximates supposing gas phase (a-aproximation)
#liquid_fzero aproximates supposing liquid phase (a-aproximation)

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



#=
Z = pv/rt 
ZRT/p

=#

function liquid_fzero(p,Pspec,t,v0 = 0.375*RGAS*t/Pspec
    ,relax=0.9*one(float(Pspec))
    ,maxevals=100
    ,atol = zero(float(Pspec))
    ,rtol = 8*eps(one(float(Pspec))))
    _Pspec,_v0,_t,_atol,_rtol,_relax = promote(Pspec,v0,t,atol,rtol,relax)
    return _liquid_fzero(p,_Pspec,_v0,_t,_relax,maxevals,_atol,_rtol)
end



#2.2V = b

function _liquid_fzero(P,Pspec::TT,v0::TT,T::TT,relax = TT(0.9),maxevals=100, atol=zero(TT), rtol=8*eps(TT)) where TT
    _1 = one(TT)
    v = v0
    count = 0
    maxcount = 10
    converged = false
    vmin = TT(Inf)
    vold=TT(NaN)
    pold = TT(NaN)
    function f0(v) 
        res = _1/P(v) - _1/Pspec
        if res < 0
           return -10*res
        else
           return res
        end
    end

    f00(v) =  P(v) - Pspec
    while !converged
        v0 = Optim.optimize(f0,eps(TT),v0,iterations=16).minimizer
        p0  = P(v0)
        if abs(vold-v0)/v0 < 0.001
            break
        end
        b1 = abs(-RGAS*T/p0 + v0)
        b2 = abs(-RGAS*T/Pspec + b1)
        vold = v0
        vmin = min(vmin,b1,b2,v0)
        v0 = vmin
        #@show v0
        if count >= maxcount
            converged = true
        end

        count +=1
        pold = p0


    end
    p0 = P(v0)
    m = ForwardDiff.derivative(P,v0)
    #m = y2-y1/(x2-x1)
    #(yx-y0) = m(xx-x0)
    vx = (Pspec - p0)/m + v0
    res =  roots_fzero(P,Pspec,T,vx)
    #@show res
    return res
end

#=
p > RT/v-b
p(v-b) > RT
v - b > RT/P
v + b> v > RT/P + b 
v> v - b> RT/P
b< v-RT/P 

=#


function v_zeros(::CubicRoots,mt::MultiPT,model,p,t,x;phase=:unspecified)
    _poly = cubic_poly(mt,model,p,t,x)
    poly = Polynomials.ImmutablePolynomial(_poly)
    sols = Polynomials.roots(poly)
    real_sols = cardan_reals(sols)
    RTp = RGAS*t/p
    vl = minimum(real_sols)*RTp
    vv = maximum(real_sols)*RTp
    if isnan(vl) & !isnan(vv)
        return (vv,vv)
    elseif isnan(vv) & !isnan(vl)
        return (vl,vl)
    else
        return (vl,vv)
    end
end

function v_zeros(::CubicRoots,mt::SinglePT,model,p,t;phase=:unspecified)
    _poly = cubic_poly(mt,model,p,t)
    poly = Polynomials.ImmutablePolynomial(_poly)
    sols = Polynomials.roots(poly)
    real_sols = cardan_reals(sols)
    RTp = RGAS*t/p
    vl = minimum(real_sols)*RTp
    vv = maximum(real_sols)*RTp
    if isnan(vl) & !isnan(vv)
        return (vv,vv)
    elseif isnan(vv) & !isnan(vl)
        return (vl,vl)
    else
        return (vl,vv)
    end
end

function v_zero(method::CubicRoots,mt::MultiPT,model,p,t,x,v0=nothing;phase=:unspecified)
    vols = v_zeros(method,mt,model,p,t,x)
    if is_liquid(phase)
        return first(vols)
    elseif is_gas(phase)
        return last(vols)
    else
        return v_zero_general(mt,model,p,t,x)
    end
end

function v_zero(method::CubicRoots,mt::SinglePT,model,p,t,v0=nothing;phase=:unspecified)
    vols = v_zeros(method,mt,model,p,t)
    if is_liquid(phase)
        return first(vols)
    elseif is_gas(phase)
        return last(vols)
    else
        return v_zero_general(mt,model,p,t)
    end
end