
#=

al cubic residual helmholtz models are from:

https://doi.org/10.6028/jres.121.011

Helmholtz Energy Transformations of Common Cubic Equations of State for Use with Pure Fluids and Mixtures


=#
function pressure_impl(mt::SingleVT,model::CubicModel,v,t)
    a,b,p = cubic_abp(mt,model,v,t)
    return p
end

function pressure_impl(mt::MultiVT,model::CubicModel,v,t,x)
    a,b,p = cubic_abp(mt,model,v,t,x)
    return p
end


"""
    cubic_mixing_rule(op, x, p,A=nothing)


returns an efficient implementation of:
`sum(x[i]*x[j]*op(i,j) for i = 1:n , j = 1:n)`

"""
function cubic_mixing_rule(op, x,A=nothing)
    if isnothing(A)
        N = length(x)
        @boundscheck checkbounds(x, N)
        @inbounds begin
            res1 = zero(eltype(x))
            for i = 1:N
                res1 += op(i,i) * x[i]^2
                for j = 1:i - 1
                    res1 += 2 * x[i] * x[j] * op(i, j)
                end
            end
        end
        return res1
    else
        N = length(x)
        @boundscheck checkbounds(x, N)
        @inbounds begin
            res1 = zero(eltype(x))
            for i = 1:N
                res1 += op(i,i) * x[i]^2
                for j = 1:i - 1
                    res1 += 2 * x[i] * x[j] * op(i, j) * A[i,j]
                end
            end
        end
        return res1
    end
end

#=
cubic polynomial:

k0 = -B((c1*c2)*B*(B+1)+A)
k1 = c1*c2*B*B-(c1+c2)*B*(B+1) + A
k2 = (c1+c2-1)*B-1
k3 = 1

on vdw:

c1=c2 = 0, c1*c2 = 0, c1+c2 = 0

k0 = -B((0)*B*(B+1)+A)
k1 = 0*B*B-(0)*B*(B+1) + A
k2 = (0-1)*B-1
k3 = 1

k0 = -AB
k1 =  A
k2 = -(B+1)
k3 = 1

on rk:

c1 = 1, c2 = 0, c1*c2 = 0,c1+c2 = 1

k0 = -B(0*B*(B+1)+A)
k1 = 0*B*B-(1)*B*(B+1) + A
k2 = (1-1)*B-1
k3 = 1

k0 = -AB
k1 = -B*(B+1) + A
k2 = -1
k3 = 1


on pr:

c1 = 1+√2, c2 = 1-√2 c1*c2 = -1,c1+c2 = 1

k0 = -B((-1)*B*(B+1)+A)
k1 = -B*B-(1)*B*(B+1) + A
k2 = (1-1)*B-1
k3 = 1


k0 = B*(B*(B+1)-A)
k1 = -2*B*B - B + A
k2 = -1
k3 = 1

k0 = B*(B*(B+1)-A)
k1 = -B*(2*B+1) + A
k2 = -1
k3 = 1

=#



