#copied and adapted from PolynomialRoots.jl
function cardan(poly::NTuple{4,T}) where {T<:AbstractFloat}
    # Cubic equation solver for complex polynomial (degree=3)
    # http://en.wikipedia.org/wiki/Cubic_function   Lagrange's method
    a1  =  1 / complex(poly[4])
    E1  = -complex(poly[3])*a1
    E2  =  complex(poly[2])*a1
    E3  = -complex(poly[1])*a1
    s0  =  E1
    E12 =  E1*E1
    A   =  2*E1*E12 - 9*E1*E2 + 27*E3 # = s1^3 + s2^3
    B   =  E12 - 3*E2                 # = s1 s2
    # quadratic equation: z^2 - Az + B^3=0  where roots are equal to s1^3 and s2^3
    Δ = sqrt(A*A - 4*B*B*B)
    if real(conj(A)*Δ)>=0 # scalar product to decide the sign yielding bigger magnitude
        s1 = exp(log(0.5 * (A + Δ)) * third)
    else
        s1 = exp(log(0.5 * (A - Δ)) * third)
    end
    if s1 == 0
        s2 = s1
    else
        s2 = B / s1
    end
    zeta1 = complex(-0.5, sqrt(T(3.0))*0.5)
    zeta2 = conj(zeta1)
    return third*(s0 + s1 + s2), third*(s0 + s1*zeta2 + s2*zeta1), third*(s0 + s1*zeta1 + s2*zeta2)
end
function volume_solver(
    mt::SinglePT,
    method::CubicRoots,
    model::HelmholtzModel,
    p,
    t,
    v0 = nothing;
    no_pts = 7)

    return 1
end