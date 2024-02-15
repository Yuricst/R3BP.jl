"""
Module related to restricted three body problem
"""
module R3BP

import Roots
using DifferentialEquations
using Printf: @printf
using LinearAlgebra
using Printf

using AstrodynamicsBase

include("misc.jl")
include("r3bp_params.jl")
include("ode/eoms.jl")

abstract type R3BPPropagator end
abstract type R3BPPropagatorSTM end
include("ode/propagator_cr3bp.jl")
include("ode/propagator_bcr4bp.jl")

include("lpo/collinear_halo.jl")
include("differentialcorrection/singleshooting.jl")
include("differentialcorrection/multipleshooting.jl")

end # module R3BP
