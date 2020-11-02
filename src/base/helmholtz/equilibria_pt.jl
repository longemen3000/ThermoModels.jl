

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


function flash_impl(mt::MultiPT,model,p,t,z0;threaded=true)
    _0 = zero(p*t*first(z0))
    _1 = one(_0)
    logϕ0 = similar(z0)
    if threaded
        res_z0 = Threads.@spawn create_logϕ!!($model,logϕ0,$p,$t,z0)
    end
    kmodel = WilsonK(model)
    k = kvalues_impl(mt,kmodel,p,t,z0)
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
    
    if threaded
        res_l = Threads.@spawn update_logϕ!!($model,logϕl,$p,$t,xil,nothing,:liquid)
        res_v = Threads.@spawn update_logϕ!!($model,logϕv,$p,$t,xiv,nothing,:vapor)
         vl,logϕl = fetch(res_l)
         vv,logϕv = fetch(res_v)
     else
         vl,logϕl = update_logϕ!!(model,logϕl,p,t,xil,nothing,:liquid)
         vv,logϕv = update_logϕ!!(model,logϕv,p,t,xiv,nothing,:vapor)
     end
     #@show logϕl,logϕv
    lnk = logϕl - logϕv
    @! k .= exp.(lnk)
  
    for i = 1:5
        res = refine_k!!(model,p,t,z0,xil,xiv,vl,vv,k,logϕl,logϕv,lnk,threaded)
        xil,xiv,vl,vv,k,logϕl,logϕv,lnk,β = res
        if !(_0 <= β <= _1)
            #@show β = β0
            break
        end
    end
    #if the process was threaded, the value is fetched here
    if !threaded
        vz0,logϕ0 = create_logϕ!!(model,logϕ0,p,t,z0)
    else
        vz0,logϕ0 = fetch(res_z0)
    end

    tpdl = flash_tpd(xil, logϕl,z0,logϕ0)
    tpdv = flash_tpd(xiv, logϕv,z0,logϕ0)
    ΔG = (1 - β) * tpdl + β * tpdv
    if all(>(0),(ΔG,tpdv,tpdl))
        stable_phase = true # more detailed stability phase analisis
        #@show stable_phase
    else
        stable_phase = false #proceed to optim
        #@show stable_phase

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
            opts = Optim.Options(
                                    iterations = 1,
                                    store_trace = true,
                                    show_trace = true)
            
        
            lower =zeros(length(z0))
            upper = Array(z0)
  
            optimize(optim_obj,nv0, Optim.LBFGS(),opts)
            
            #nv and nl results can be complete garbage,
            #but xil and xiv are good.
            @show vle_res.xil
            @show vle_res.xiv
            β = vfrac_from_fracs(z0,vle_res.xil,vle_res.xiv)
            @show β
            return vtx_states(vle_res)
    end
    #return original phase, for now an
    #todo: identify phase
    return state(t=t,v=vz0,xn = z0)
    #"return error("phase is not stable")
    #println(stable_phase)
    #==))
    interres =  Dict{Symbol,Any}(:xl =>xil
    ,:xv => xiv
    ,:vl=> vl
    ,:vv =>vv
    ,:k=>k
    ,:logϕl=> logϕl
    ,:logϕv=> logϕv
    ,:lnk=>lnk
    ,:stable_phase=>stable_phase)
    ==#
    vle_res = VLEResults(
    model = model
    ,z0 = z0 #initial composition
    ,vl = vl #liquid mol volume
    ,vv = vv #gas mol volume
    ,tv = t #temperature
    ,pv = p #pressure
    ,tl = t #temperature
    ,pl = p #pressure
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
    #println(nv0)
    #println(fg!(1,nothing,nv0))

    optim_obj = OnceDifferentiable(only_fg!(fg!), nv0)
    opts = Optim.Options()
    opts = Optim.Options(
                             iterations = 10,
                             store_trace = true,
                             show_trace = true)
    

    optimize(optim_obj, nv0, Optim.BFGS(),opts)
    
    return vtx_states(vle_res)
end



function refine_k!!(model,p,t,z0,xil,xiv,vl,vv,k,logϕl,logϕv,lnk,threaded)  
    if threaded
       res_l = Threads.@spawn update_logϕ!!(model,logϕl,$p,$t,xil,$vl,:liquid)
       res_v = Threads.@spawn update_logϕ!!(model,logϕv,$p,$t,xiv,$vv,:vapor)
        vl,logϕl = fetch(res_l)
        vv,logϕv = fetch(res_v)
    else
        vl,logϕl = update_logϕ!!(model,logϕl,p,t,xil,vl,:liquid)
        vv,logϕv = update_logϕ!!(model,logϕv,p,t,xiv,vv,:vapor)
    end
    @! lnk .= logϕl .- logϕv
    @! k .= exp.(lnk)
  
    β = flash_vfrac(RR(),k,z0) 
    xil = flash_liquid!!(xil,k,z0,β)
    xiv = flash_vapor!!(xiv,k,z0,β)
    xil = normalizefrac!!(xil)
    xiv = normalizefrac!!(xiv)
    
    



    return xil,xiv,vl,vv,k,logϕl,logϕv,lnk,β
    
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


#create value of logϕ for stream composition. selects volume with the lowest gibbs energy
function create_logϕ!!(model,logϕ,p,t,x,v=nothing)
    x = normalizefrac!!(x)
    _vtx = QuickStates.vtx()      
    _ptx = QuickStates.ptx()   
    v =  v_zero(_ptx,model,p,t,x,v)
    @! logϕ .= fugacity_coeff_impl(_vtx,model,v,t,x)
    return v,logϕ
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

function vtx_states(obj::VLEResults)
    stᵥ = state(mol_v=obj.vv,t=obj.tv,xn=obj.xiv,phase=:gas)
    stₗ = state(mol_v=obj.vl,t=obj.tl,xn=obj.xil,phase=:liquid)
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
        
  
        
  
        obj.nl = obj.z0 - obj.nv

        
        obj.nl = remove_minimums!!(obj.nl)
        obj.nv = remove_minimums!!(obj.nv)
        if !iszero(sum(obj.nl))
            obj.xil = normalizefrac!!(obj.nl)
        end
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
    return flash_impl(mt,model,p,t,z0,threaded=thr)
end

#v = xiv*β
#l = z0 - v

#dfugv = I/β - 1/vt + dlnfugv/vt
#dfugl = I/(1-β) - 1/lt + dlnfugl/lt