

#from PC-SAFT paper https://doi.org/10.1021/ie0003887
function fugacity_coeff_impl(mt::MultiVT,model::HelmholtzModel,v,t,x)
    Z = compresibility_factor_impl(mt,model,v,t,x)
    ar = z -> αR_impl(mt,model,one(v)/v,t,z)
    dadx = ForwardDiff.gradient(ar,x)
    sumdadx = dot(x,dadx)
    lnZ = log(Z)
    c0 = ar(x) + (Z-1) - lnZ - sumdadx
    return map(xi->xi+c0,dadx)
end

function fugacity_coeff_impl(mt::SingleVT,model::HelmholtzModel,v,t)
    rho = one(v)/v
    ar(_rho) = αR_impl(mt,model,_rho,t)
    (ar,dadrho) = f_dfdx(ar,rho)
    p = rho*(one(rho)+rho*dadrho)*RGAS*t
    Z = p*v/(RGAS*t)
    return ar + (Z-1) - log(Z)
end

function flash_impl(mt::MultiPT,model,p,t,z0)
    function v_solver(p0,t0,x0,v0=nothing,pts=7) 
        return volume_solver(QuickStates.ptx(),volume_solver_type(model),model,p0,t0,x0,v0,no_pts=pts)
    end
    _vtx = QuickStates.vtx()      
    kmodel = WilsonK(model)
    k = kvalues_impl(mt,kmodel,p,t)
    (vl,vg) = v_solver(p,t,z0)
    g_0 = flash_eval(k,z0,0.0)
    g_1 = flash_eval(k,z0,1.0)
    if g_0 <= 0  #bubble point assumption

        xil = z0
        xiv = flash_vapor(k,z0,0.0)
    elseif g_1 >= 0 #dew point assumption
        xil = flash_liquid(k,z0,1.0)
        xiv = z0
    else #two phase assumption
        β = flash_vfrac(RR(),k,z0)
        xil = flash_liquid(k,z,β)
        xiv = flash_vapor(k,z,β)
    end
    println(xil)
    println(xiv)

    xil = normalizefrac!!(xil)
    xiv = normalizefrac!!(xiv)
    _vl = Threads.@spawn v_solver(p,t,xil)
    _vv = Threads.@spawn v_solver(p,t,xiv)
    vl = first(fetch(_vl))
    vv = last(fetch(_vv))
    @show logϕl = fugacity_coeff_impl(_vtx,model,vl,t,xil)
    @show logϕv = fugacity_coeff_impl(_vtx,model,vv,t,xiv)
    lnk = logϕl - logϕv
    k .= exp.(lnk)
    β = flash_vfrac(RR(),k,z0)
end

