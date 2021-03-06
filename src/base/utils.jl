const RGAS = 8.31446261815324

unicode_subscript(a::Int) = join('₀'+d for d in reverse(digits(a)))

function molar_to_weight(value,x,mw)
    return value/(mw*x)*1000.0
end

# mixing rule:
#= =
an operation on a vector of properties that returns a single sumber.
is a weighted operation on p, by the weights x.
# 

= =#
"""
    mixing_rule(op, x, p)

returns an efficient implementation of:
`sum(x[i]*x[j]*op(p[i],p[j]) for i = 1:n , j = 1:n)`
where `op(p[i],p[j]) == op(p[j],p[i])`

""" 
function mixing_rule(op, x, p,A=nothing)
    if isnothing(A)
        N = length(x)
        @boundscheck checkbounds(p, N)
        @inbounds begin
            res1 = zero(eltype(x))
            for i = 1:N
                res1 += p[i] * x[i]^2
                for j = 1:i - 1
                    res1 += 2 * x[i] * x[j] * op(p[i], p[j])
                end
            end
        end
        return res1
    else
        return mixing_rulea(op,x,p,A)
    end
end
# 
# Abstract mixing rule, generates a mixing rule,based on 
# an operation, so mij = xi*xj*op(pi,pj)*Aij
# example: mixing_rule(geometric_mean_rule,x,Tc,1.-K)
# 

function mixing_rulea(op, x, p, A)
    N = length(x)
    checkbounds(A, N, N)
    @boundscheck checkbounds(p, N)
    @inbounds begin
        res1 = zero(eltype(x))
        for i = 1:N
            res1 += p[i] * x[i]^2
            for j = 1:i - 1
                res1 += 2 * x[i] * x[j] * op(p[i], p[j]) * A[i, j]
            end
        end
    end
    return res1
end
# 
# Abstract asymetric mixing rule 
# Adds a Asymetric matrix A_asym, and a op_sim(xi,xj,Aij)
# the mayor example is the GERG2008 equation, where
# op_asym(xi,xj,Aij) = (xi+xj)/(Aij^2*xi + xj)
# 
"""
    mixing_rule_asymetric(op, op_asym, x, p, A, A_asym)

returns an efficient implementation of:
` sum(A[i,j] * x[i] * x[j] * op(p[i],p[j]) * op_asym(x[i],x[j],A_asym[i,j])) for i = 1:n , j = 1:n)`
where `op(p[i],p[j]) == op(p[j],p[i])` , op_asym doesn't follow this symmetry.

""" 
function mixing_rule_asymetric(op, op_asym, x, p, A, A_asym)
    N = length(x)
    checkbounds(A, N, N)
    checkbounds(A_asym, N, N)
    @boundscheck checkbounds(p, N)
    @inbounds begin
        res1 = zero(eltype(x))

        for i = 1:N
            x[i] != 0 && begin
                res1 += p[i] * x[i]^2
                for j = 1:i - 1
                    res1 +=
                        2 *
                        x[i] *
                        x[j] *
                        op(p[i], p[j]) *
                        A[i, j] *
                        op_asym(x[i], x[j], A_asym[i, j])
                end
            end
        end
    end
    return res1
end



"""
    mixing_matrix!(A, op, p)

returns a matrix of size nxn , where A[i,j] = op(p[i],p[j]), modifying A inplace

""" 
function mixing_matrix!(A, op, p)
    N = length(p)
    @boundscheck size(A) == (N, N)
    @inbounds begin
        res1 = zero(eltype(p))
        for i = 1:N
            A[i, i] = p[i]
            for j = 1:i - 1
                A[i, j] = op(p[i], p[j])
                A[j, i] = op(p[i], p[j])
            end
        end
    end
    return A

end







"""
    mixing_matrix!(op, p)

returns a matrix of size nxn , where A[i,j] = op(p[i],p[j])
""" 
function mixing_matrix(op, p)
    N = length(p)
    A = Matrix{eltype(p)}(undef, N, N)
    return mixing_matrix!(A, op, p)
end

#a_in_b(t1, t2) = all(in(t2), t1)
#tuple_comparison(t1, t2) = ((length(t1) == length(t2)) && a_in_b(t1, t2) && a_in_b(t2, t1))


#returns f(x)/(∂f(x)/∂x), for Newton solvers, using ForwardDiff+StaticArrays+DiffResults for one single AD pass

"""
    autonewton(f,x)

returns f/(df/fx) evaluated in `x`, using `ForwardDiff.jl`, `DiffResults.jl` and `StaticArrays.jl` to calculate f and dfdx in one pass
"""
function autonewton(f,x::T) where T
    _f(z) = f(only(z))
    x_vec =   SVector(x)
    ∂result = DiffResults.GradientResult(x_vec)  
    _∂f =  ForwardDiff.gradient!(∂result, _f,x_vec)
    fx =  DiffResults.value(_∂f)
    ∂f∂x = only(DiffResults.gradient(_∂f))
    return fx/∂f∂x
end

function f_dfdx(f,x::T) where T
    _f(z) = f(only(z))
    x_vec =   SVector(x)
    ∂result = DiffResults.GradientResult(x_vec)  
    _∂f =  ForwardDiff.gradient!(∂result, _f,x_vec)
    fx =  DiffResults.value(_∂f)
    ∂f∂x = only(DiffResults.gradient(_∂f))
    return (fx,∂f∂x)
end
"""
    normalizefrac!!(x)

mutates x, so that sum(x) == 1. also removes negative values and replaces Inf with big numbers

"""
function normalizefrac!!(x)
    x0 = zero(eltype(x))
    function f(xi) 
        if xi < x0
         return -x0 
        elseif isinf(xi)
            return inv(eps(typeof(xi))) 
        else
            return xi
        end
    end
    xn = one(x0)/sum(f,x)
    @! x .= f.(x) .* xn
    return x
end


function nan_to_num(xi)
    if x == typemax(xi)
        return prevfloat(xi)
    elseif isnan(xi)
        return zero(xi)
    else
        return xi
    end
end

function remove_minimums!!(x,lo=zero(eltype(x)))
    return @! x .= max.(x,lo)
end

function remove_maximums!!(x,hi=inv(eps(eltype(x))))
    return @! x .= min.(x,hi)
end

function cardan_reals(res)
    nan = real(zero(first(res))/zero(first(res)))
    a1, a2, a3 = res
    isreal1 = abs(imag(a1) < 4*eps(imag(a1)))
    isreal2 = abs(imag(a2) < 4*eps(imag(a2)))
    isreal3 = abs(imag(a3) < 4*eps(imag(a3)))
    r1 = real(a1)
    r2 = real(a2)
    r3 = real(a3)
    if isreal1 & isreal2 & isreal3
        return (r1,r2,r3)
    elseif !isreal1 & !isreal2
        return (r3,r3,r3)
    elseif  !isreal3 & !isreal2
        return (r1,r1,r1)
    elseif  !isreal3 & !isreal1
        return (r2,r2,r2)
    elseif !isreal1
        return (r3,r2,r2)
    elseif !isreal2
        return (r1,r3,r3)
    elseif !isreal3
        return (r1,r3,r3)
    else
        return (nan,nan,nan)
    end
end


function is_liquid(x::Symbol)::Bool
    x in (:liquid,:l,:L,:LIQUID)
end

function is_gas(x::Symbol)::Bool
    x in (:gas,:g,:G,:GAS,:vapor,:VAPOR,:V,:v)
end

