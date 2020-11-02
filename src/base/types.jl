

abstract type AbstractVolumeSolver end
struct VolumeBisection <: AbstractVolumeSolver 
    pts::Int64
end
VolumeBisection() = VolumeBisection(21)

struct CubicRoots <: AbstractVolumeSolver end
struct SUVA <: AbstractVolumeSolver end


abstract type HelmholtzModel <: ThermoState.ThermoModel end
abstract type CubicModel <: HelmholtzModel end

abstract type SaturationModel <: ThermoState.ThermoModel end
abstract type GasModel <: ThermoState.ThermoModel end
abstract type SatLiquidModel <: ThermoState.ThermoModel end

#catch all to build models that doesnt need automatically
#to stop erroring



#abstract types corresponding 
#to certain points on the thermodynamic params space

primitive type MyEnum <: Base.Enums.Enum{Int32} 32 end 
function MyEnum(x::Integer)
    1<=x<=3 || Base.Enums.enum_argument_error(MyEnum,x)
    return Base.bitcast(MyEnum, convert(Int32, x))
end

Base.Enums.namemap(x::Type{MyEnum}) = Dict{Int32,Symbol}(1=>:a,2=>:b,3=>:c)

Base.typemin(x::Type{MyEnum}) = MyEnum(1)
Base.typemax(x::Type{MyEnum}) = MyEnum(3)

const a = MyEnum(1)
const b = MyEnum(2)
const c = MyEnum(3)
Base.instances(x::Type{MyEnum}) = (a,b,c)
