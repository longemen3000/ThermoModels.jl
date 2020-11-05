

#from PC-SAFT paper https://doi.org/10.1021/ie0003887
function fugacity_coeff_impl(mt::MultiVT,model::HelmholtzModel,v,t,x)
    Z = compressibility_factor_impl(mt,model,v,t,x)
    ar = z -> αR_impl(mt,model,one(v)/v,t,z)
    dadx = ForwardDiff.gradient(ar,x)
    sumdadx = dot(x,dadx)
    lnZ = log(Z)
    c0 = ar(x) + (Z-1) - lnZ - sumdadx
    return map(xi->xi+c0,dadx)
end


function flash_impl(mt::MultiPT,model,p,t,z0;threaded=true,moles=one(eltype(z0)))
    _0 = zero(p*t*first(z0))
    _1 = one(_0)
    logϕ0 = similar(z0)
    if threaded
        res_z0 = Threads.@spawn update_logϕ!!(model,logϕ0,$p,$t,z0,nothing,:undefined)
    end
    kmodel = WilsonK(model)
    k0 = kvalues_impl(mt,kmodel,p,t,z0)
    k = copy(k0)
    #@show k
    g_0 = flash_eval(k,z0,_0)
    g_1 = flash_eval(k,z0,_1)
    if g_0 <= 0  #bubble point assumption
        #println("bubble point")
        β = _0
        β0 = _0
        xil = copy(z0)
        xiv = flash_vapor(k,z0,_0)
    elseif g_1 >= 0 #dew point assumption
        #println("dew point")
        β = _1
        β0 = _1
        xil = flash_liquid(k,z0,_1)
        xiv = copy(z0)
    else #two phase assumption
        #println("flash point")
        β = flash_vfrac(RR(),k,z0)
        xil = flash_liquid(k,z0,β)
        xiv = flash_vapor(k,z0,β)
    end
    logϕl = similar(xil)
    logϕv = similar(xiv)
    xil = normalizefrac!!(xil)
    xiv = normalizefrac!!(xiv)
    
    vl,logϕl,vv,logϕv = update_logϕ_vl!!(model,p,t,z0,xil,xiv,logϕl,logϕv,nothing,nothing,threaded)

    lnk = logϕl - logϕv
    @! k .= exp.(lnk)
    stable_phase = false
    for i = 1:5
        res = refine_k!!(model,p,t,z0,xil,xiv,vl,vv,k,logϕl,logϕv,lnk,threaded)
        xil,xiv,vl,vv,k,logϕl,logϕv,lnk,β = res
        if !(_0 <= β <= _1)
            stable_phase = true
            break
        end
    end
    #if the process was threaded, the value is fetched here
    if !threaded
        vz0,logϕ0 = update_logϕ!!(model,logϕ0,p,t,z0,nothing,:undefined)
    else
        vz0,logϕ0 = fetch(res_z0)
    end

    if stable_phase == false
        tpdl = flash_tpd(xil, logϕl,z0,logϕ0)
        tpdv = flash_tpd(xiv, logϕv,z0,logϕ0)
        ΔG = (1 - β) * tpdl + β * tpdv
        if all(>(0),(ΔG,tpdv,tpdl))
            stable_phase = true # more detailed stability phase analisis
        end
    end

    if stable_phase == true
        mul_v = copy(xiv)
        mul_l = copy(xiv)
        xil,xiv,vl,vv = trial_compositions!!(model,p,t,z0,xil,xiv,vl,vv,logϕl,logϕv,logϕ0,mul_v,mul_l,true;k0=k0)  
        count = 1
        max_count = 10
        while count < max_count
            xil,xiv,vl,vv,tpdl,tpdv = trial_compositions!!(model,p,t,z0,xil,xiv,vl,vv,logϕl,logϕv,logϕ0,mul_v,mul_l,true)  
            converged,stable_phase = trial_convergence(z0,xil,xiv,mul_l,mul_v,tpdl,tpdv,stable_phase)
            if converged == true
                break
            else
                count += 1
            end
        end
        if stable_phase == false
            for i in 1:3
                xil,xiv,vl,vv,tpdl,tpdv = trial_compositions!!(model,p,t,z0,xil,xiv,vl,vv,logϕl,logϕv,logϕ0,mul_v,mul_l,true)  
            end
        end
    end

    if stable_phase == false
        vle_res = VLEResults(
            model = model
            ,vl = vl #liquid mol volume
            ,vv = vv #gas mol volume
            ,tv = t #temperature
            ,pv = p #pressure
            ,tl = t #temperature
            ,pl = p #pressure
            ,z0 = z0 #initial composition
            ,lnϕl = logϕl #efective log fugacity coef - liquid
            ,lnϕv = logϕv #efective log fugacity coef - gas
            ,nl = copy(xil) #moles of liquid
            ,nv = copy(xiv) #moles of gas
            ,xil = xil #mol fraction of liquid
            ,xiv = xiv) #mol fraction of gas
        
            fg! = generate_ptflash_objective(vle_res,threaded=threaded)
            nv0 = zeros(length(xiv))
            nv0 .+= xiv
            nv0 .*= β
        optim_obj = OnceDifferentiable(only_fg!(fg!), nv0)
        optimize(optim_obj,nv0, Optim.LBFGS())
        return vtx_states(vle_res,moles)
    end
  
    return state(t=t,v=vz0,xn = z0)
end


function trial_compositions!!(model,p,t,z0,xil,xiv,vl,vv,logϕl,logϕv,logϕ0,mul_v,mul_l,threaded=true;k0=nothing)  
    if k0 !== nothing
        @! xiv .= z0 .* k0
        @! xil .= z0 ./ k0
        xil = normalizefrac!!(xil)
        xiv = normalizefrac!!(xiv) 
        vl = nothing
        vv = nothing
        tpdv = convert(typeof(p),Inf)
        tpdl = convert(typeof(p),Inf)
    else
        vl,logϕl,vv,logϕv = update_logϕ_vl!!(model,p,t,z0,xil,xiv,logϕl,logϕv,nothing,nothing,threaded)
        mul_v = exp.(logϕ0 .- logϕv) 
        mul_l = exp.(logϕ0 .- logϕl)
        xiv = z0 .* mul_v
        xil = z0 .* mul_l
        xil = normalizefrac!!(xil)
        xiv = normalizefrac!!(xiv) 
        tpdl = flash_tpd(xil, logϕl,z0,logϕ0)
        tpdv = flash_tpd(xiv, logϕv,z0,logϕ0)
       
    end
    return xil,xiv,vl,vv,tpdl,tpdv
end

function trial_convergence(z0,xil,xiv,mul_l,mul_v,tpdl,tpdv,stable_phase)
    new_stable_phase = stable_phase
    @! mul_l .= mul_l .- one(eltype(mul_l)) 
    @! mul_v .= mul_v .- one(eltype(mul_v)) 
    #this checks for:
    # a: the composition aproaches z0
    # b: variation between iterations is very small

    if norm(mul_v) ≈ 0.0
        converged = true
     elseif norm(mul_l) ≈ 0.0
         converged = true
     elseif xiv ≈ z0
        converged = true
    elseif xil ≈ z0
        converged = true
    elseif tpdl < 0
        converged = true
        new_stable_phase = false
    elseif tpdv < 0
        converged = true
        new_stable_phase = false
    else
        converged = false
    end
    return converged,new_stable_phase
end

function refine_k!!(model,p,t,z0,xil,xiv,vl,vv,k,logϕl,logϕv,lnk,threaded)  
    vl,logϕl,vv,logϕv = update_logϕ_vl!!(model,p,t,z0,xil,xiv,logϕl,logϕv,vl,vv,threaded)
    @! lnk .= logϕl .- logϕv
    @! k .= exp.(lnk)
  
    β = flash_vfrac(RR(),k,z0) 
    xil = flash_liquid!!(xil,k,z0,β)
    xiv = flash_vapor!!(xiv,k,z0,β)
    xil = normalizefrac!!(xil)
    xiv = normalizefrac!!(xiv)

    return xil,xiv,vl,vv,k,logϕl,logϕv,lnk,β
    
end

#update both liquid and gas fugacity coefficients, considering threading
function update_logϕ_vl!!(model,p,t,z0,xil,xiv,logϕl,logϕv,vl=nothing,vv=nothing,threaded=true)  
    if threaded
        res_l = Threads.@spawn update_logϕ!!(model,logϕl,$p,$t,xil,$vl,:liquid)
        res_v = Threads.@spawn update_logϕ!!(model,logϕv,$p,$t,xiv,$vv,:vapor)
         vl,logϕl = fetch(res_l)
         vv,logϕv = fetch(res_v)
     else
         vl,logϕl = update_logϕ!!(model,logϕl,p,t,xil,vl,:liquid)
         vv,logϕv = update_logϕ!!(model,logϕv,p,t,xiv,vv,:vapor)
     end
     return vl,logϕl,vv,logϕv

end
#with values of p,t,x; returns a new logϕ,v
function update_logϕ!!(model,logϕ,p,t,x,v=nothing,phase=:liquid)
    x = normalizefrac!!(x)
    _vtx = QuickStates.vtx()    
    _ptx = QuickStates.ptx()   
    _v =  v_zero(_ptx,model,p,t,x,v;phase=phase)
    @! logϕ = fugacity_coeff_impl(_vtx,model,_v,t,x)
    return _v,logϕ
end


function flash_tpd(xi, lnϕ,x0,lnϕ0)
    f(xiᵢ,lnϕᵢ,x0ᵢ,lnϕ0ᵢ) = xiᵢ*(log(xiᵢ) + lnϕᵢ-log(x0ᵢ)- lnϕ0ᵢ)
    return mapreduce(f,+,xi, lnϕ,x0,lnϕ0)
end


Base.@kwdef mutable struct VLEResults{MODEL,S,V}
    model::MODEL
    vl::S #liquid mol volume
    vv::S #gas mol volume
    tv::S #temperature, gas
    pv::S #pressure, gas
    tl::S #temperature,liquid
    pl::S #pressure,liquid
    z0::V #initial composition
    lnϕl::V #efective log fugacity coef - liquid
    lnϕv::V#efective log fugacity coef - gas
    nl::V #moles of liquid
    nv::V #moles of gas
    xil::V #mol fraction of liquid
    xiv::V #mol fraction of gas
end

function vtx_states(obj::VLEResults,moles=1.0)
    β = sum(obj.nv)
    nv = β*moles
    nl = (1-β)*moles
    stᵥ = state(mol_v=obj.vv,t=obj.tv,xn=obj.xiv,moles=nv,phase=:gas)
    stₗ = state(mol_v=obj.vl,t=obj.tl,xn=obj.xil,moles=nl,phase=:liquid)
    return (stₗ,stᵥ)
end

function vtn_states(obj::VLEResults)
    stᵥ = state(mol_v=obj.vv,t=obj.tv,n=obj.nv,phase=:gas)
    stₗ = state(mol_v=obj.vl,t=obj.tl,n=obj.nl,phase=:liquid)
    return (stₗ,stᵥ)
end

function _flash_obj(xlᵢ,xvᵢ,lnϕlᵢ,lnϕvᵢ,β)
    _0 = zero(β)
    l = ifelse(xlᵢ>_0,(1-β)*xlᵢ*(log(xlᵢ)+lnϕlᵢ),_0)
    v = ifelse(xvᵢ>_0,β*xvᵢ*(log(xvᵢ)+lnϕvᵢ),_0)
    return l+v
end

function generate_ptflash_objective(obj::VLEResults;threaded = true)
    #TODO: implement threading
    function flash_objective!(F, G, nv)
        t = obj.tv
        p = obj.pv
        β = sum(nv)
        obj.nv = nv
        obj.nl = abs.(obj.z0 - obj.nv)
        #corr = minimum(obj.nl)
        #@show obj.nl
        #obj.nl = obj.z0 - obj.nv - corr
        #@show obj.nl
        obj.nl = remove_minimums!!(obj.nl)
        obj.nv = remove_minimums!!(obj.nv)
        obj.xil = normalizefrac!!(obj.nl)
        obj.xiv = normalizefrac!!(obj.nv)

        #@show obj.xil
        obj.vl ,obj.lnϕl = update_logϕ!!(obj.model,obj.lnϕl,p,t,obj.xil,obj.vl,:liquid)
        obj.vv ,obj.lnϕv = update_logϕ!!(obj.model,obj.lnϕv,p,t,obj.xiv,obj.vv,:gas)
    
        if !(G === nothing)
            G .= zero(eltype(G))
            G .+= obj.lnϕv .+ log.(obj.xiv)   .- obj.lnϕl .- log.(obj.xil) 
        end
        if !(F === nothing)
            fobjective = (xlᵢ,xvᵢ,lnϕlᵢ,lnϕvᵢ)->_flash_obj(xlᵢ,xvᵢ,lnϕlᵢ,lnϕvᵢ,β) 
            F = mapreduce(fobjective,+,obj.xil,obj.xiv,obj.lnϕl,obj.lnϕv)
            return F
        end
    end
    return flash_objective!
end

function equilibria(mt::MultiPT,model::HelmholtzModel,st::ThermodynamicState)
    t = temperature(FromState(),st)
    opts = options(FromState(),st)
    thr = get(opts,:threaded,true)
    p = pressure(FromState(),st)
    z0 = mol_fraction(FromState(),st,nothing,molecular_weight(model))
    _moles = moles(FromState(),st,u"mol",molecular_weight(model))
    return flash_impl(mt,model,p,t,z0,threaded=thr,moles=_moles)
end

#v = xiv*β
#l = z0 - v

#dfugv = I/β - 1/vt + dlnfugv/vt
#dfugl = I/(1-β) - 1/lt + dlnfugl/lt