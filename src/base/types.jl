
const SINGLE = SingleComponent
const MULTI = MaterialCompounds{MOLAR,FRACTION}

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

