
abstract type KValueModel <: ThermoModel end

function kvalues(model::KValueModel,st::ThermodynamicState)
    return kvalues(state_type(st),model,st)
end



abstract type RachfordRiceSolver end
struct RachfordRice <: RachfordRiceSolver end
struct MultiPhaseRachfordRice <: RachfordRiceSolver end
struct PolynomialRachfordRice <: RachfordRiceSolver end
struct LeiboviciNichita  <: RachfordRiceSolver end
struct LiJohnsAhmadi <: RachfordRiceSolver end

struct AnalyticalRachfordRiceSolution{N} <: RachfordRiceSolver end


const RR = RachfordRice
const MultiRR = MultiPhaseRachfordRice
const PolynomialRR = PolynomialRachfordRice
const LN2 = LeiboviciNichita
const LJA = LiJohnsAhmadi



function vf0(K,z)
    N = length(K)
    Kmin = minimum(K)
    Kmax,i_kmax = findmax(K)
    z_of_Kmax = z[i_kmax]
    Kminm1 = one(Kmin)-Kmin
    V_over_F_min = ((Kmax-Kmin)*z_of_Kmax - Kminm1)/(Kminm1*(Kmax - one(Kmin)))
    V_over_F_max = one(Kmin)/Kminm1
    return (V_over_F_min,V_over_F_max)
end




function rr_find_zero(f0,(vmin,vmax),method::Roots.AbstractBracketing)
    return Roots.find_zero(f0,(vmin,vmax),method)
end

function rr_find_zero(f0,(vmin,vmax),method::Roots.AbstractNonBracketing)
    β0 = (vfmin+vfmax)/2
    return Roots.find_zero(f0,β0,method)
end


function rr_find_zero(f0,(vmin,vmax),method::Roots.Newton)
    β0 = (vfmin+vfmax)/2
    AD_f0(x) = autonewton(f0,x)
    return Roots.find_zero(f0,β0,method)
end

function flash_vfrac(model::RR,K,Z,roots_method=Roots.Order0())
    Kmin0 = minimum(K)
    Kmax0 = maximum(K)
    if (Kmin0 > 1.0) | (Kmax0 < 1.0)
        return ArgumentError("no positive composition solution with given K-values")
    end
    vfmin,vfmax = vf0(K,Z)

    function rr(k,z,β)
        km1 = k-one(k)
        return z*km1/(1+β*km1)
    end
    f0(β) = mapreduce((k,z)->rr(k,z,β),+,K,Z)
    βres = rr_find_zero(f0,(vfmin,vfmax),method)
    return res
end

function flash_vfrac(model::LeiboviciNichita,K,Z,roots_method=Roots.Newton())
    Kmin0 = minimum(K)
    Kmax0 = maximum(K)
    if (Kmin0 > 1.0) | (Kmax0 < 1.0)
        return ArgumentError("no positive composition solution with given K-values")
    end
    vfmin,vfmax = vf0(K,Z)
    ck = vfmin
    ckm1 = vfmax
    _λ(y) = ck + (ckm1-ck)/(one(y)+exp(-y))
    _y(λ) = log((ckm1-ck)/(λ-ck)-one(λ))
    function rr(k,z,y)
        c = 1/(1-k)
        return z/_y((y)-c)
    end

    f0(β) = mapreduce((k,z)->rr(k,z,β),+,K,Z)
    yres = rr_find_zero(f0,(_y(vfmin),_y(vfmax)),roots_method)
    return _λ(yres)
end

function flash_vfrac(model::PolynomialRachfordRice,K,Z,method=nothing)
    n = length(Z)
    if 2 <= n <= 4
        return flash_vfrac(AnalyticalRachfordRiceSolution{n}(),K,Z)
    else
        return DomainError("analytical solutions are only available for 2,3,or 4 compounds")
    end
end


function rr0(k,z)
    km1 = k-one(k)
    return z*km1
end 

function rr1(k,z)
    k1 = one(k)
    return z*(k1-k1/k)
end

function rry(k,z,y)
    k1 = one(k)
    km1 = k-k1
    return z*km1/(k1+β*km1)
end         

function flash_eval(K,Z,β)  
    if iszero(β)
        return mapreduce((k,z)->rr0(k,z),+,K,Z)
    elseif isone(β)
        return mapreduce((k,z)->rr1(k,z),+,K,Z)
    else
        return mapreduce((k,z)->rry(k,z,β),+,K,Z)
    end
end



struct WilsonK{T} <: ThermoModel 
    pc::T
    tc::T
    ω::T
end

struct MollerupK{T} <: ThermoModel
    pc::T
    tc::T
    ω::T
end



function kvalues(mt::MultiPT,model::WilsonK,st::ThermodynamicState)
    p = pressure(FromState(),st)
    t = temperature(FromState(),st)
    return kvalues_impl(mt,model,p,t)
end

function kvalues_impl(mt::MultiPT,model::WilsonK,p,t)
    pc = pressure(model,CriticalPoint())
    tc = temperature(model,CriticalPoint())
    ω = acentric_factor(model)
    return exp.(log.(pc./p).+5.373 .*(1.0 .+ ω).*(1.0 .-tc./t))
end

function kvalues(mt::MultiPT,model::MollerupK,st::ThermodynamicState)
    p = pressure(FromState(),st)
    t = temperature(FromState(),st)
    return kvalues_impl(mt,model,p,t)
end

function kvalues_impl(mt::MultiPT,model::MollerupK,p,t)
    pc = pressure(model,ThermoState.CriticalPoint())
    tc = temperature(model,ThermoState.CriticalPoint())
    ω = acentric_factor(model)
    return  (pc ./ p) .* exp.(5.42 .* (1.0 .- (tc ./ t)))
end

