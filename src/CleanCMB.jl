module CleanCMB
using LinearAlgebra
export harmonic_ilc_weights, harmonic_ilc_clean_cl
export harmonic_cilc_weights, harmonic_cilc_clean_cl
include("harmonic_ilc.jl")
include("harmonic_cilc.jl")
end
