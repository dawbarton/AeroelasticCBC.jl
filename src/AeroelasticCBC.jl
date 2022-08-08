module AeroelasticCBC

using StaticArrays
using ComponentArrays
using BifurcationKit
using Setfield
using LinearAlgebra
using PGFPlotsX

include("model.jl")
include("bifurcation.jl")
include("plots.jl")

end
