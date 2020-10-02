"""
    Region4

    Returns the vapour-liquid phase boundary (saturation line).
    Valid from triple point to critical point
    273.15K ≤ T ≤ 647.096K
    InputType is either :T or :P to indicate that InputValue is temperature [K]
    or pressure [MPa]. The complimentary value is returned.
"""
function Region4 end

const Region4_tup = (;n = [0.116_705_214_527_67E4,
    -0.724_213_167_032_06E6,
    -0.170_738_469_400_92E2,
    0.120_208_247_024_70E5,
    -0.323_255_503_223_33E7,
    0.149_151_086_135_30E2,
    -0.482_326_573_615_91E4,
    0.405_113_405_420_57E6,
    -0.238_555_575_678_49,
    0.650_175_348_447_98E3]
)

function pressure_impl(mt::SingleSatT,model::IF97Region{:r4},t)
    n = Region4_tup.n
    Θ = t + n[9]/(t - n[10])
    A =      Θ^2 + n[1]*Θ + n[2]
    B = n[3]*Θ^2 + n[4]*Θ + n[5]
    C = n[6]*Θ^2 + n[7]*Θ + n[8]
    P = (2C / (-B + √(B^2 - 4*A*C)))^4
    return convert_unit(u"MPa",u"Pa",P)
end

function temperature_impl(mt::SingleSatP,model::IF97Region{:r4},p)
    n = Region4_tup.n
    P = convert_unit(u"Pa",u"MPa",p)
    β = P^0.25
    E =      β^2 + n[3]*β + n[6]
    F = n[1]*β^2 + n[4]*β + n[7]
    G = n[2]*β^2 + n[5]*β + n[8]
    D = 2G / (-F - √(F^2 - 4*E*G))
    T = (n[10]+D-√((n[10]+D)^2 - 4(n[9]+n[10]*D)))/2.0
    return T
end

pressure_impl(mt::SingleΦT,model::IF97Region{:r4},t) = pressure_impl(QuickStates.sat_p(),model,t)
temperature_impl(mt::SingleΦP,model::IF97Region{:r4},p) = temperature_impl(QuickStates.sat_t(),model,p)

