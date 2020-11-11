include("vdw.jl")
include("rk.jl")
include("rks.jl")
include("pr.jl")
const VdW = VanDerWaals
const RK = RedlichKwong
const PR = PengRobinson
const SRK = SoaveRedlichKwong

export VdW
export RK
export SRK
export PR
