module AeroelasticCBC

using StaticArrays: @SVector, @SMatrix
using ComponentArrays: ComponentArrays, ComponentArray
using BifurcationKit: BifurcationKit, BifurcationProblem, ContinuationPar, NewtonPar,
                      PeriodicOrbitTrapProblem, PALC, continuation, ContinuousEvent, SaveAtEvent
using Setfield: @lens, @set, @set!
using LinearAlgebra: norm
using OrdinaryDiffEq: ODEProblem, solve, Tsit5
using PGFPlotsX: PGFPlotsX, @pgf, Axis, PlotInc

const OUTPATH = joinpath(@__DIR__, "output")

include("model.jl")
include("bifurcation.jl")
include("tartaruga/Tartaruga.jl")
include("plotting.jl")

end # module

# enlarge_x_limits = false
