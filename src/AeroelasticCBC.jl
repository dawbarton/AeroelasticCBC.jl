module AeroelasticCBC

using StaticArrays: @SVector, @SMatrix
using ComponentArrays: ComponentArrays, ComponentArray
using BifurcationKit: BifurcationKit, BifurcationProblem, ContinuationPar, NewtonPar,
                      PeriodicOrbitTrapProblem, PALC, continuation, ContinuousEvent
using Setfield: @lens, @set, @set!
using LinearAlgebra: norm
using OrdinaryDiffEq: ODEProblem, solve, Tsit5
using PGFPlotsX: PGFPlotsX, @pgf

include("model.jl")
include("bifurcation.jl")
include("plots.jl")

end
