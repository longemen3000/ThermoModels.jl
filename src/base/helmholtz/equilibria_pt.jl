

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
    p = float(p)
    t = float(t)
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
        initial_data = (model,p,t,z0)
        xil,xiv,vl,vv = trial_compositions!!(initial_data,xil,xiv,vl,vv,logϕl,logϕv,logϕ0,mul_v,mul_l,true;k0=k0)  
        count = 1
        max_count = 10
        while count < max_count
            xil,xiv,vl,vv,tpdl,tpdv = trial_compositions!!(initial_data,xil,xiv,vl,vv,logϕl,logϕv,logϕ0,mul_v,mul_l,true)  
            converged,stable_phase = trial_convergence(z0,xil,xiv,mul_l,mul_v,tpdl,tpdv,stable_phase)
            if converged == true
                break
            else
                count += 1
            end
        end
        if stable_phase == false
            for i in 1:3
                xil,xiv,vl,vv,tpdl,tpdv = trial_compositions!!(initial_data,xil,xiv,vl,vv,logϕl,logϕv,logϕ0,mul_v,mul_l,true)  
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
            #@show vle_res.xiv
            #@show vle_res.xil
            #@show vle_res.nv
            #@show vle_res.nl
            fg! = generate_ptflash_objective(vle_res,threaded=threaded)
            nv0 = zeros(length(xiv))
            nv0 .+= xiv
            nv0 .*= β
            opts = Optim.Options(g_tol = 1e-12,
                iterations = 10,
                store_trace = true,
                show_trace = true)
            optim_obj = OnceDifferentiable(only_fg!(fg!), nv0)
         try
            optimize(optim_obj,nv0, Optim.BFGS(linesearch = Optim.LineSearches.MoreThuente()))
         catch err
            @show vle_res
         finally
            return vtx_states(vle_res,moles)
         end
    end
  
    return state(t=t,v=vz0,xn = z0)
end


function trial_compositions!!(initial_data,xil,xiv,vl,vv,logϕl,logϕv,logϕ0,mul_v,mul_l,threaded=true;k0=nothing)  
    model,p,t,z0 = initial_data
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
        @! xiv += z0 .* mul_v - xiv
        @! xil += z0 .* mul_l - xil
        xil = normalizefrac!!(xil)
        xiv = normalizefrac!!(xiv) 
        tpdl = flash_tpd(xil, logϕl,z0,logϕ0)
        tpdv = flash_tpd(xiv, logϕv,z0,logϕ0)
        @show tpdl
        @show tpdv
       
    end
    return xil,xiv,vl,vv,tpdl,tpdv
end

function trial_convergence(z0,xil,xiv,mul_l,mul_v,tpdl,tpdv,stable_phase)
    new_stable_phase = stable_phase
    _ones = ones(length(mul_l))
  
    #this checks for:
    # a: the composition aproaches z0
    # b: variation between iterations is very small

    if mul_v ≈ _ones
        converged = true
     elseif mul_l ≈ _ones
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
    #@show xil,xiv,vl,vv
    @! lnk .= logϕl .- logϕv
    @! k .= exp.(lnk)
   k = remove_maximums!!(k)
    β = flash_vfrac(RR(),k,z0) 
    #@show β
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
     #@show logϕl
     #@show logϕv
     #@! lnk .= logϕl .- logϕv
     #@! k .= exp.(lnk)
     #@show k
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
    xv = obj.xiv
    xl = obj.xil
    z0 = obj.z0
    β = (z0[1] - xl[1])/(xv[1]-xl[1])
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
        #@show nv
        obj.nl =obj.z0 - obj.nv
        #corr = minimum(obj.nl)
        #@show obj.nv
        #obj.nl = obj.z0 - obj.nv - corr
        #@show obj.nl
        obj.nl = remove_minimums!!(obj.nl)
        obj.nv = remove_minimums!!(obj.nv)
        obj.xil = normalizefrac!!(obj.nl)
        obj.xiv = normalizefrac!!(obj.nv)
        
        
        #@show obj.xil
        obj.vl ,obj.lnϕl = update_logϕ!!(obj.model,obj.lnϕl,p,t,obj.xil,obj.vl,:liquid)
        obj.vv ,obj.lnϕv = update_logϕ!!(obj.model,obj.lnϕv,p,t,obj.xiv,obj.vv,:gas)
        #@show obj.lnϕl
        #@show obj.lnϕv
        if !(G === nothing)
            #@show G
            G .= zero(eltype(G))
            for i in 1:length(G)
                if obj.xiv[i] <= eps(obj.xiv[i])
                    G[i] += obj.lnϕv[i]  - log(obj.xil[i]) 
                elseif obj.xil[i] <= eps(obj.xil[i])
                    G[i] += log(obj.xiv[i])   - obj.lnϕl[i]
                else
                    G[i] += obj.lnϕv[i] + log(obj.xiv[i])   - obj.lnϕl[i] - log(obj.xil[i]) 
                #G .+= obj.lnϕv .+ log.(obj.xiv)   .- obj.lnϕl .- log.(obj.xil) 
                end
            #@show G
            end
        end
        if !(F === nothing)
            fobjective = (xlᵢ,xvᵢ,lnϕlᵢ,lnϕvᵢ)->_flash_obj(xlᵢ,xvᵢ,lnϕlᵢ,lnϕvᵢ,β) 
            F = mapreduce(fobjective,+,obj.xil,obj.xiv,obj.lnϕl,obj.lnϕv)
            #@show F
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

#=

water-o2 (GERG2008)
logϕl = [-0.7705794332079066, 14.528869934800472]
logϕv = [-0.009360929108700455, 0.0031516545028221]
k = [0.46709692059495184, 2.034413807942005e6]

logϕl = [-0.7705794332079066, 14.528869934800472]
logϕv = [-0.009360929108700452, 0.0031516545028221037]
k = [0.46709692059495184, 2.034413807942005e6]

logϕl = [-0.7705816360739632, 14.539131948004911]
logϕv = [-0.007992505036874167, 0.0017762880364641798]
k = [0.4664571435209862, 2.0582273486877447e6]

logϕl = [-0.770581636063973, 14.539132013140865]
logϕv = [-0.007980392574881982, 0.0017656849290196598]
k = [0.4664514936154415, 2.0582493064751942e6]

logϕl = [-0.7705816360640014, 14.53913201314133]
logϕv = [-0.007980285727848792, 0.0017655915176444346]
k = [0.4664514437764727, 2.0582494987400582e6]

logϕl = [-0.7705816360640014, 14.53913201314133]
logϕv = [-0.007980284785338597, 0.0017655906936611884]
k = [0.4664514433368374, 2.05824950043602e6]

water-o2 (VdW)
logϕl = [2.297372915993078, 7.406207150097631]
logϕv = [-0.005105694210269807, -6.038786149844145e-5]
k = [9.998935228783932, 1646.2702241435495]

logϕl = [2.297372915993078, 7.406207150097631]
logϕv = [-0.005105694210269807, -6.038786149844145e-5]
k = [9.998935228783932, 1646.2702241435495]

logϕl = [2.297372728601342, 7.406704158470847]
logϕv = [-0.003828630603820233, -0.0005244515961122089]
k = [9.986172231292917, 1647.8531690781929]

logϕl = [2.297411298402489, 7.389040025797274]
logϕv = [-0.005202637502164627, 6.585518575033793e-5]
k = [10.000288433265684, 1618.0454227403882]

logϕl = [2.297412692566457, 7.388698667718173]
logϕv = [-0.005203259181312256, 6.678678750513239e-5]
k = [10.000308592298598, 1617.4916772651154]
=#

