module CleanCMB
using LinearAlgebra
export ilc_weights, ilc_clean_cij
export cilc_weights
export milca_weights
include("ilc.jl")
include("cilc.jl")
include("milca.jl")
end
