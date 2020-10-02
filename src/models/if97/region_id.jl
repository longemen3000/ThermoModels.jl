function B23(InputType::Pressure, InputValue)
    n = B23_n
    return n[4] + √((InputValue - n[5])/n[3])
end

function B23(InputType::Temperature, InputValue)
    n = B23_n
    return n[1] + n[2] * InputValue + n[3] * InputValue^2
end

const B2bc_n = [0.905_842_785_147_23E3,
                -0.679_557_863_992_41,
                0.128_090_027_301_36E-3,
                0.265_265_719_084_28E4,
                0.452_575_789_059_48E1]
"""
    B2bc

    Returns the boundary between regions 2b and 2c, approx s=5.85kJ/kgK
    InputType is either :h or :P to indicate that InputValue is enthalpy [kJ/kg]
    or pressure [MPa]. The complimentary value is returned.
    Valid from saturation line at 554.485K & 6.54670MPa to 1019.32K & 100MPa
"""
function B2bc(InputType::Enthalpy, InputValue)
    n = B2bc_n
    return n[1] + n[2] * InputValue + n[3] * InputValue^2
end

function B2bc(InputType::Pressure, InputValue)
    n = B2bc_n
    return n[1] + n[2] * InputValue + n[3] * InputValue^2
end


function temperature_impl(mt::SingleSatP,model::WaterIF97,p)
P = normalize_units(p)
Pc = pressure(model,CriticalPoint())
Pmin = 611.213 #different than triple point pressure
    if Pmin <= P <= Pc
        return temperature_impl(mt,IF97Region{:r4}(),P)
    else
        return DomainError(p, "Pressure not between triple and critical points")
    end
end


#temperature_impl(mt::SingleSatP,model::WaterIF97,p)
function pressure_impl(mt::SingleSatT,model::WaterIF97,t)
    T = normalize_units(t)
    Tc = temperature(model,CriticalPoint())
    T3= temperature(model,TriplePoint())
    if T3 <= T <= Tc
        return pressure_impl(mt,IF97Region{:r4}(),T)
    else
        return DomainError(t, "temperature not between triple and critical points")
    end 
end

pressure_impl(mt::SingleΦT,model::WaterIF97,t) = pressure_impl(QuickStates.sat_p(),model,t)
temperature_impl(mt::SingleΦP,model::WaterIF97,p) = temperature_impl(QuickStates.sat_t(),model,p)

function pressure(model::WaterIF97,st::ThermodynamicState,unit = u"Pa")
    return pressure(state_type(st),model,st,unit)
end

function pressure(mt::SingleSatT,model::WaterIF97,st::ThermodynamicState,unit)
    t = temperature(FromState(),st)
    res = pressure_impl(mt,model,t)
    return convert_unit(u"Pa",unit,res)
end

function temperature(model::WaterIF97,st::ThermodynamicState,unit = u"K")
    return pressure(state_type(st),model,st,unit)
end

function temperature(mt::SingleSatP,model::WaterIF97,st::ThermodynamicState,unit)
    t = pressure(FromState(),st)
    res = temperature_impl(mt,model,t)
    return convert_unit(u"K",unit,res)
end



function region_id(mt::SinglePT,model::WaterIF97,p, t)::Symbol
    function PSat(Tx)
        p = pressure_impl(QuickStates.sat_t(),model,Tx) 
        return convert_unit(u"Pa",u"MPa",p)
    end
    TSat(Px) =  temperature_impl(QuickStates.sat_t(),model,Px) 
    P = convert_unit(u"Pa",u"MPa",p)
    T = normalize_units(t)
    #=
    Region 1:
        273.15K ≤ T ≤ 623.15K  Psat(T) ≤ P ≤ 100MPa

    Region 2:
    273.15K ≤ T ≤ 623.15K  0 ≤ P ≤ Psat(T)
    623.15K <  T ≤ 863.15K  0 <  P ≤ P(T) from B23-model
    863.15K <  T ≤ 1073.15K 0 <  P ≤ 100MPa

    Region 2meta:
    From the saturated-vapour line to the 5% equlibrium moisture line, a.k.a the
    practical Wilson line:
        %equilibrium moisture = [h - h_liq(P)] / [h_vap(P) - h_liq(P)]
    for pressures from the triple point (273.16K, 611.657MPa) up to 10MPa.
    These equations only used on demand.

    Region 3:
        623.15K ≤ T ≤ T(P) from B23-model P(T) from B23-model ≤ P ≤ 100MPa

    Region 4:
    Returns the vapour-liquid phase boundary (saturation line).
    Valid from triple point to critical point
    273.15K ≤ T ≤ 647.096K
    Since this is a single floating point number, Region 4 is never returned and only used on demand.

    Region 5:emperature in [K]
    =#

    #Check region 5 first:
    if 1073.15 ≤ T ≤ 2273.15 && 0 ≤ P ≤ 50
        return :r5
    elseif T < 273.15 || T > 1073.15 || P < 0 || P > 100
        throw(DomainError((p, t), "Pressure/Temperature outside valid ranges.")) #Outside the regions where equations are available.
    elseif 273.15 ≤ T ≤ 623.15
        if Psat(T) ≤ P ≤ 100
            return :r1
        else
            return :r2
        end
    elseif 623.15 <  T ≤ 863.15
        if 0 <  P ≤ B23(Temperature(), T)
            return :r2
        else
            return :r3
        end
    else
        return :r2
    end
end
function region_id(mt::SinglePH,model::WaterIF97,p, _h)::Symbol
    
    function PSat(Tx)
        p = pressure_impl(QuickStates.sat_t(),model,Tx) 
        return convert_unit(u"Pa",u"MPa",p)
    end
    TSat(Px) =  temperature_impl(QuickStates.sat_t(),model,Px) 
    P = convert_unit(u"Pa",u"MPa",p)
    h = convert_unit(u"kJ/g",u"kJ/kg",_h)
    #=
        Region 1:
        273.15K ≤ T ≤ 623.15K  Psat(T) ≤ P ≤ 100MPa

        Region 2:
        273.15K ≤ T ≤ 623.15K  0 ≤ P ≤ Psat(T)
        623.15K <  T ≤ 863.15K  0 <  P ≤ P(T) from B23-model
        863.15K <  T ≤ 1073.15K 0 <  P ≤ 100MPa
        For reverse:
            2a: P ≤ 4Mpa
            2b: s ≥ 5.85 kJ/kgK (or use function B2bc if P,h specificed)
            2c: s ≥ 5.85 kJ/kgK (or use function B2bc if P,h specificed)
    =#

    # Check overall region first
    if P > 100.0
        throw(DomainError(P, "Pressure not in valid ranges."))
    end

    T = temperature_impl(mt,IF97Region{:r1}(),p, _h)
    if 273.15 ≤ T ≤ 623.15
        # could be Region 1
        if P ≥ Psat(T)
            return :r1 #, Region1(:SpecificH, P, T) # Return forward h for consistency check
        end
    end

    if P ≤ 4.0
        # could be Region 2a
        T = temperature_impl(mt,IF97Region{:r2a}(),p, _h)
        if Tsat(P) ≤ T ≤ 1073.15
            return :r2a #, Region2(:SpecificH, P, T) # Return forward h for consistency check
        else
            throw(DomainError((P, T), "Pressure/Temperature outside valid ranges.")) # Only other region with backwards mdoels for P ≤ 4.0 is Region 1, which is already eliminated
        end
    else
        htest = B2bc(Pressure(), P)
        if h < htest
            # could be Region 2c
            T = temperature_impl(mt,IF97Region{:r2c}(),p, _h)
            if P ≤ Psat(623.15)
                if Tsat(P) ≤ T ≤ 1073.15
                    return :r2c #, Region2(:SpecificH, P, T) # Return forward h for consistency check
                end
            elseif B23(Pressure(), P) ≤ T ≤ 1073.15
                return :r2c #, Region2(:SpecificH, P, T) # Return forward h for consistency check
            else
                throw(DomainError((P, T), "Pressure/Temperature outside valid ranges."))
            end
        else
            # could be Region 2b
            T = temperature_impl(mt,IF97Region{:r2b}(),p, _h)
            if P ≤ Psat(623.15)
                if Tsat(P) ≤ T ≤ 1073.15
                    return :r2b #, Region2(:SpecificH, P, T) # Return forward h for consistency check
                end
            elseif B23(Pressure(), P) ≤ T ≤ 1073.15
                return :r2b #, Region2(:SpecificH, P, T) # Return forward h for consistency check
            else
                throw(DomainError((P, T), "Pressure/Temperature outside valid ranges."))
            end
        end
    end
    throw(DomainError((P, T), "Pressure/Temperature outside valid ranges."))
end

function region_id(mt::SinglePS,model::WaterIF97,p, _s)::Symbol
    function PSat(Tx)
        p = pressure_impl(QuickStates.sat_t(),model,Tx) 
        return convert_unit(u"Pa",u"MPa",p)
    end
    TSat(Px) =  temperature_impl(QuickStates.sat_t(),model,Px) 

    P = convert_unit(u"Pa",u"MPa",p)
    s = convert_unit(u"J/(kg*K)",u"kJ/(kg*K)",_s)
    #=
        Region 1:
        273.15K ≤ T ≤ 623.15K  Psat(T) ≤ P ≤ 100MPa

        Region 2:
        273.15K ≤ T ≤ 623.15K  0 ≤ P ≤ Psat(T)
        623.15K <  T ≤ 863.15K  0 <  P ≤ P(T) from B23-model
        863.15K <  T ≤ 1073.15K 0 <  P ≤ 100MPa
        For reverse:
            2a: P ≤ 4Mpa
            2c: s ≥ 5.85 kJ/kgK
            2b: s ≥ 5.85 kJ/kgK
    =#

    # Check overall region first
    if P > 100
        throw(DomainError((P, s), "Pressure/entropy not in valid ranges."))
    end

    T = temperature_impl(mt,IF97Region{:r1}(),p, _s)

    if 273.15 ≤ T ≤ 623.15
        # could be Region 1
        if P ≥ Psat(T)
            return :r1 #, Region1(:SpecificS, P, T) # Return forward s for consistency check
        end
    end

    if P ≤ 4.0
        # could be Region 2a
        T = temperature_impl(mt,IF97Region{:r2a}(),p, _s)
        if Tsat(P) ≤ T ≤ 1073.15
            return :r2a #, Region2(:SpecificS, P, T) # Return forward s for consistency check
        else
            throw(DomainError((P, s), "Pressure/entropy not in valid ranges.")) # Only other region with backwards mdoels for P ≤ 4.0 is Region 1, which is already eliminated
        end
    else
        if s < 5.85
            # could be Region 2c
            T = temperature_impl(mt,IF97Region{:r2c}(),p, _s)
            if P ≤ Psat(623.15)
                if Tsat(P) ≤ T ≤ 1073.15
                    return :r2 #c, Region2(:SpecificS, P, T) # Return forward s for consistency check
                end
            elseif B23(Pressure(),P) ≤ T ≤ 1073.15
                return :r2c #, Region2(:SpecificS, P, T) # Return forward s for consistency check
            else
                throw(DomainError((P, s), "Pressure/entropy not in valid ranges."))
            end
        else
            # could be Region 2b
            T = temperature_impl(mt,IF97Region{:r2b}(),p, _s)
            if P ≤ Psat(623.15)
                if Tsat(P) ≤ T ≤ 1073.15
                    return :r2b #, Region2(:SpecificS, P, T) # Return forward s for consistency check
                end
            elseif B23(Pressure(), P) ≤ T ≤ 1073.15
                return :r2b #, Region2(:SpecificS, P, T) # Return forward s for consistency check
            else
                throw(DomainError((P, s), "Pressure/entropy not in valid ranges."))

            end
        end
    end
    throw(DomainError((P, s), "Pressure/entropy not in valid ranges."))
end
