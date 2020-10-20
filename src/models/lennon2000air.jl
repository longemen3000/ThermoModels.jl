struct Lennnon2000Air <: ThermoState.ThermoModel end
const TAU_MAX_EXP_87 = 0.4207493606569795
const lemmon2000_air_R = 8.314510

const lemmon2000_air_T_reducing = 132.6312
const lemmon2000_air_P_reducing = 3.78502E6
const lemmon2000_air_rho_reducing = 10447.7
const lemmon2000_air_rho_reducing_inv = 1.0/lemmon2000_air_rho_reducing

const lemmon2000_air_MW = 28.9586
const lemmon2000_air_P_max = 2000E6
const lemmon2000_air_T_max = 2000.

molecular_weight(::Lennnon2000Air) = lemmon2000_air_MW

function _f0(::Lennnon2000Air, delta,tau)
    tau_inv = one(tau)/tau
    A0 =  (-0.00019536342*tau*sqrt(tau) + 17.275266575*tau + tau_inv*(tau_inv*(6.057194e-8*tau_inv 
            - 0.0000210274769) - 0.000158860716) + log(delta) + 2.490888032*log(tau)
    
            # These two logs both fail for tau < 1e-18, can be truncated but should not be necessary.
            + 0.791309509*log(1.0 - exp(-25.36365*tau)) + 0.212236768*log(1.0 - exp(-16.90741*tau)) 
            - 13.841928076)
    
            if tau < TAU_MAX_EXP_87
        A0 -= 0.197938904*log(exp(87.31279*tau) + (2.0/3.0))
    else
        A0 -= 17.282597957782162*tau # 17.282... = 87.31279*0.197938904
    return A0
    end
end

function _fr(::Lennnon2000Air,delta,tau)

    delta2 = delta*delta
    delta3 = delta*delta2
    delta4 = delta2*delta2
    delta5 = delta*delta4
    delta6 = delta2*delta4
    taurt2 = sqrt(tau)
    taurt4 = sqrt(taurt2)
    tau2 = tau*tau
    tau3 = tau*tau2
    tau6 = tau3*tau3
    tau12 = tau6*tau6
    tau_100 = tau^0.01
    tau2_100 = tau_100*tau_100
    tau4_100 = tau2_100*tau2_100
    tau5_100 = tau_100*tau4_100
    tau10_100 = tau5_100*tau5_100
    tau15_100 = tau5_100*tau10_100
    tau8_100 = tau4_100*tau4_100
    tau16_100 = tau8_100*tau8_100
    tau20_100 = tau4_100*tau16_100
    tau32_100 = tau16_100*tau16_100
    tau33_100 = tau_100*tau32_100
    tau64_100 = tau32_100*tau32_100
    tau80_100 = tau16_100*tau64_100
    tau40_100 = tau20_100*tau20_100
    tau97_100 = tau33_100*tau64_100
    tau45_100 = tau5_100*tau40_100
    tau90_100 = tau45_100*tau45_100
    tau160_100 = tau80_100*tau80_100
    tau320_100 = tau160_100*tau160_100
    x0 = exp(-delta)
    x1 = exp(-delta2)
    x2 = tau3*exp(-delta3)
    return (-0.101365037911999994*delta*tau160_100*x0 + 0.713116392079000017*delta*tau33_100
            - 0.146629609712999986*delta*tau40_100*x1*tau320_100 
            - 1.61824192067000006*delta*tau4_100*tau97_100 + 0.0148287891978000005*delta*taurt2*x2 
            + 0.118160747228999996*delta + 0.0714140178971000017*delta2 + 0.134211176704000013*delta3*tau15_100 
            - 0.031605587982100003*delta3*tau6*x1 - 0.17381369096999999*delta3*tau80_100*x0 
            - 0.00938782884667000057*delta3*x2*tau12 - 0.0865421396646000041*delta3 - 0.042053322884200002*delta4*tau20_100 
            + 0.0349008431981999989*delta4*tau2_100*tau33_100 + 0.0112626704218000001*delta4
            - 0.0472103183731000034*delta5*tau15_100*x0*tau80_100 + 0.000233594806141999996*delta5*tau3*taurt4*x1*delta6 
            - 0.0122523554252999996*delta6*tau*taurt4*x0 + 0.000164957183186000006*delta6*tau45_100*tau90_100)
end


function αR_impl(mt::SingleVT,::Lennnon2000Air, _rho, T)
    R = lemmon2000_air_R
    delta = rho*lemmon2000_air_rho_reducing_inv
    tau = lemmon2000_air_T_reducing/T
    return _fr(model, delta, tau) 
end

function α0_impl(mt::SingleVT,::Lennnon2000Air, _rho, T)
    R = lemmon2000_air_R
    delta = rho*lemmon2000_air_rho_reducing_inv
    tau = lemmon2000_air_T_reducing/T
    return _f0(model, delta, tau) 
end


function mol_helmholtzR_impl(mt::SingleVT,::Lennnon2000Air, v, t)
    rho = 1.0e-3 / v
    R = lemmon2000_air_R
    delta = lemmon2000_air_rho_reducing_inv/v
    tau = lemmon2000_air_T_reducing/T
    return R*t*_fr(model, delta, tau) 
end

function mol_helmholtz0_impl(mt::SingleVT,::Lennnon2000Air, v, t)
    rho = 1.0e-3 / v
    R = lemmon2000_air_R
    delta = lemmon2000_air_rho_reducing_inv/v
    tau = lemmon2000_air_T_reducing/T
    return R*t*_f0(model, delta, tau) 
end

function mol_helmholtz_impl(mt::SingleVT,::Lennnon2000Air, v, t)
    rho = 1.0e-3 / v
    R = lemmon2000_air_R
    delta = lemmon2000_air_rho_reducing_inv/v
    tau = lemmon2000_air_T_reducing/T
    return R*t*(_f0(model, delta, tau)+_fr(model, delta, tau))
end

