const IF97_R =  0.461526
const IF97_Tc = 647.096      #K          Critical temperature of water
const IF97_Pc = 22.064       #Pa         Critical pressure of water
const IF97_ρc = 322.         #kg/m3      Critical density of water
const IF97_T3 = 273.16       #K          Triple point temperature of water
const IF97_P3 = 611.657E-6   #MPa        Triple point pressure of water
const IF97_Mr = 18.01528     #kg/kmol    Molecular weight of water

struct WaterIF97 <: ThermoModel end

temperature(model::WaterIF97,st::CriticalPoint,unit=u"K")  = convert_unit(u"K",unit,647.096)
pressure(model::WaterIF97,st::CriticalPoint,unit=u"Pa")  = convert_unit(u"Pa",unit,2.2064e7)
mol_volume(model::WaterIF97,st::CriticalPoint,unit=u"m^3/mol")  = convert_unit(u"m^3/mol",unit,5.594803726708074e-5)
mol_density(model::WaterIF97,st::CriticalPoint,unit=u"mol/m^3")  = convert_unit(u"mol/m^3",unit,17873.72799560906)
mass_density(model::WaterIF97,st::CriticalPoint,unit=u"kg/m^3")  = convert_unit(u"kg/m^3",unit,322.0)
mass_volume(model::WaterIF97,st::CriticalPoint,unit=u"m^3/kg")  = convert_unit(u"m^3/kg",unit,0.003105590062111801)
temperature(model::WaterIF97,st::TriplePoint,unit=u"K")  = convert_unit(u"K",unit,273.16)
temperature(model::WaterIF97,st::NormalBoilingPoint,unit=u"K")  = convert_unit(u"K",unit,373.15)
pressure(model::WaterIF97,st::TriplePoint,unit=u"Pa")  = convert_unit(u"Pa",unit,611.657)
molecular_weight(model::WaterIF97) = 18.01528

const B23_n = [0.348_051_856_289_69E3,
                -0.116_718_598_799_75E1,
                0.101_929_700_393_26E-2,
                0.572_544_598_627_46E3,
                0.139_188_397_788_70E2]
"""
    B23(::Union{ThermoState.Types.Pressure,ThermoState.Types.Temperature}, InputValue)

Returns the boundary between regions 2 and 3.
InputType is either :T or :P to indicate that InputValue is temperature [K] or pressure [MPa]. The complimentary value is returned.
"""

struct IF97Region{T} <: ThermoModel end



include("region1.jl")
include("region2.jl")
include("region3.jl")
include("region4.jl")
include("region5.jl")
include("region_id.jl")



#SinglePT
    for op in [:helmholtz, :gibbs, :internal_energy, :enthalpy,:cp,:cv,:volume,:entropy]
        mol_op_impl = Symbol(:mol_,op,:_impl)
        mass_op_impl = Symbol(:mass_,op,:_impl)
        total_op_impl = Symbol(:total_,op,:_impl)
        mol_op = Symbol(:mol_,op)
        mass_op = Symbol(:mass_,op)
        total_op = Symbol(:total_,op)
        if op == :volume
            _unit = u"m^3/kg"
            mol_unit = u"m^3/mol"
            total_unit = u"m^3"

        elseif op in (:cv,:cp,:entropy)
            _unit = u"J/(kg*K)"
            mol_unit = u"J/(mol*K)"
            total_unit = u"J/(K)"
        else
            _unit = u"J/(kg)"
            mol_unit = u"J/(mol)"
            total_unit = u"J"
        end
     
        @eval begin
            function $mass_op_impl(mt::SinglePT,model::WaterIF97,p,t)
                id = region_id(mt,model,p,t)
                return $mass_op_impl(mt,IF97Region{id}(),val1,val2)
            end


            function $mass_op(model::WaterIF97,st::ThermodynamicState,unit=$_unit)
                return $mass_op(mt,model,st,unit)
            end

            function $mass_op(mt::SinglePT,model::WaterIF97,st::ThermodynamicState,unit)
                p = pressure(FromState(),st)
                t = temperature(FromState(),st)
                return convert_unit($_unit,unit,$mass_op_impl(mt,model,p,t))
            end

            function $mol_op(model::WaterIF97,st::ThermodynamicState,unit=$mol_unit)
                prod = molar_mass(FromState(),st,u"kg/mol",molecular_weight(model))
                res =  $mass_op(mt,model,st,unit)*prod
                return convert_unit($mol_unit,unit,res)
            end 
        end

        if !(op in (:cv,:cp))
            
            @eval begin
                function $total_op(model::WaterIF97,st::ThermodynamicState,unit=$total_unit)
                    prod = mass(FromState(),st,u"kg",molecular_weight(model))
                    res =  $mass_op(mt,model,st,unit)*prod
                    return convert_unit($total_unit,unit,res)
                end
            end
        end

    if op != :enthalpy
        @eval begin
            function $mass_op_impl(mt::SinglePH,model::WaterIF97,p,t)
                id = region_id(mt,model,p,t)
                return $mass_op_impl(mt,IF97Region{id}(),val1,val2)
            end

            function $mass_op(mt::SinglePH,model::WaterIF97,st::ThermodynamicState,unit)
                p = pressure(FromState(),st)
                t = temperature(FromState(),st)
                return convert_unit($_unit,unit,$mass_op_impl(mt,model,p,t))
            end
        end

            
    end
end

#==
for op in (:helmholtz, :gibbs, :internal_energy, :enthalpy)
    mol_op_impl = Symbol(:mol_,op,:_impl)
    mass_op_impl = Symbol(:mass_,op,:_impl)
    total_op_impl = Symbol(:total_,op,:_impl)
    mol_op = Symbol(:mol_,op)
    mass_op = Symbol(:mass_,op)
    total_op = Symbol(:total_,op)
    for spec in [:SinglePT,:SinglePH,:SinglePS]
    @eval begin 
        function $mass_op_impl(mt,model::HelmholtzModel,p,t)
            id = region_id(mt,model,p,t)
            return $mass_op_impl()
        end 

        function ThermoState.$mol_op(model::HelmholtzModel,st::ThermodynamicState,unit=u"J/mol")
            return $mol_op(state_type(st),model,st,unit)
        end    

        function $mol_op(mt::SingleVT,model::HelmholtzModel,st::ThermodynamicState,unit = u"J/mol")
            mw = molecular_weight(model)
            v = mol_volume(FromState(),st,u"m^3/mol",mw)
            t = temperature(FromState(),st)
            val = $mol_op_impl(mt,model,v,t)
            return convert_unit(u"J/mol",unit,val)
        end

    end
end

==#
