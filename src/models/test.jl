
using LinearAlgebra, ComponentArrays,Optim
Optim.ldiv!(out::ComponentVector, P::Optim.InverseDiagonal, A::ComponentVector) = copyto!(out, A .* P.diag)

function test()
k_unr_out = [ 0.,  3.,  1.,  5.,  2.,  4.,  5.,  0.,  0.,  2.,  0.,  4.,  5.,
         6.,  1., 22.,  0.,  1., 19.,  0.,  0.,  7.,  0.,  4.,  0.,  6.]

k_unr_in = [ 7.,  4.,  9.,  3.,  2.,  3.,  4., 11.,  4.,  4.,  5.,  4.,  3.,
         3.,  3.,  0.,  2.,  6.,  0.,  3.,  0.,  1.,  8.,  2.,  3.,  3.]

k_recip = [ 1.,  8.,  2.,  8.,  3.,  7.,  6.,  1.,  1.,  7.,  4.,  5.,  4.,
         2.,  9.,  3.,  2.,  6.,  6.,  2., 25.,  6.,  4.,  3.,  3.,  6.]

N = length(k_unr_out)

v0 = ComponentArray(x=repeat([0.1], N), y=repeat([0.1], N), z=repeat([0.1], N))

function neg_llhood(v, k_unr_out, k_unr_in, k_recip)
    res1 = (k,vx) -> mapreduce(+,(x,y)->x*log(y),k,vx)
    llhood = res1(k_unr_out,v.x)
    llhood += res1(k_unr_in,v.y)
    llhood += res1(k_recip,v.z)
    
        
    xy = v.x .* v.y'
    zz = v.z .* v.z'

    Q = log.(1 .+ xy + xy' + zz)
    
    llhood -= sum(tril(Q, -1))
    
    -llhood
end

f = v -> neg_llhood(v, k_unr_out, k_unr_in, k_recip)


lower = ComponentArray(x=repeat([floatmin()], N), y=repeat([floatmin()], N), z=repeat([floatmin()], N))
upper = ComponentArray(x=repeat([Inf], N), y=repeat([Inf], N), z=repeat([Inf], N))
inner_optimizer = LBFGS()
opts = Optim.Options(outer_iterations = 100,iterations=20)
results = optimize(f, lower, upper, v0, Fminbox(inner_optimizer),autodiff=:forward)

return results
end