module CleanCMB
using LinearAlgebra
export ilc_weights, ilc_clean_cij
export cilc_weights, cilc_clean_cij
include("ilc.jl")
include("cilc.jl")
end
