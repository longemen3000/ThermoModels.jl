function calcz(b0, c0)
    z0 = (1.0 + sqrt(1.0 + 4*b0)) / 2.0
    f(y) = -c0 + y*(-b0 + y*(-1.0 + y))
    return Roots.find_zero(f,z0)   
end
