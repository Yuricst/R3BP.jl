"""
Module related to restricted three body problem
"""
module R3BP

import Roots
using DifferentialEquations
using Printf: @printf

using AstrodynamicsBase

include("r3bp_params.jl")
include("ode/eoms.jl")

abstract type R3BPPropagator end
include("ode/propagator_cr3bp.jl")
include("ode/propagator_bcr4bp.jl")

end # module R3BP
