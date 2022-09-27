module AeroelasticCBC

using StaticArrays: @SVector, @SMatrix
using ComponentArrays: ComponentArrays, ComponentArray
using BifurcationKit: BifurcationKit, BifurcationProblem, ContinuationPar, NewtonPar,
                      PeriodicOrbitTrapProblem, PALC, continuation, ContinuousEvent
using Setfield: @lens, @set, @set!
using LinearAlgebra: norm
using OrdinaryDiffEq: ODEProblem, solve, Tsit5
using PGFPlotsX: PGFPlotsX, @pgf

push!(PGFPlotsX.CUSTOM_PREAMBLE, "\\usepackage{unicode-math}")
push!(PGFPlotsX.CUSTOM_PREAMBLE, "\\setmainfont{TeX Gyre Heros}")
push!(PGFPlotsX.CUSTOM_PREAMBLE, "\\setmathfont{TeX Gyre Pagella Math}")
push!(PGFPlotsX.CLASS_OPTIONS, "12pt")

include("model.jl")
include("bifurcation.jl")
include("tartaruga/Tartaruga.jl")


end # module

# enlarge_x_limits = false
