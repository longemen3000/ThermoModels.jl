include("vdw.jl")
include("rk.jl")
include("rks.jl")
include("pr.jl")
const VdW = VanDerWaals
const RK = RedlichKwong
const RKS = RedlichKwongSoave
const PR = PengRobinson

export VdW
export RK
export RKS
export PR
