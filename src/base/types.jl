

abstract type AbstractVolumeSolver end
struct VolumeBisection <: AbstractVolumeSolver end
struct CubicRoots <: AbstractVolumeSolver end


abstract type HelmholtzModel <: ThermoState.ThermoModel end
abstract type CubicModel <: HelmholtzModel end

abstract type SaturationModel <: ThermoState.ThermoModel end
abstract type GasModel <: ThermoState.ThermoModel end
abstract type SatLiquidModel <: ThermoState.ThermoModel end



#abstract types corresponding 
#to certain points on the thermodynamic params space


