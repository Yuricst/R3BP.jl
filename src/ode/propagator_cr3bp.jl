"""
Propagator classes with R3BP
"""


"""
Propagator for CR3BP
"""
mutable struct CR3BPPropagator <: R3BPPropagator
    problem
    method
    reltol::Float64
    abstol::Float64

    function CR3BPPropagator(
        method,
        mu::Real;
        reltol::Float64=1e-12,
        abstol::Float64=1e-12
    )
        # construct ODE problem
        problem = ODEProblem(
            rhs_cr3bp_sv!,
            [1.0, 0.0, 0.0, 0.5, 1.0, 0.0],  # placeholder for u0
            [0.0, 1.0],                      # placeholder for tspan
            (mu,),
        )

        # initialize
        new(
            problem,
            method,
            reltol,
            abstol,
        )
    end
end


"""
Propagator for CR3BP with STM
"""
mutable struct CR3BPPropagatorSTM <: R3BPPropagator
    problem
    method
    reltol::Float64
    abstol::Float64

    function CR3BPPropagator(
        method,
        mu::Real;
        reltol::Float64=1e-12,
        abstol::Float64=1e-12
    )
        # construct ODE problem
        problem = ODEProblem(
            rhs_cr3bp_svstm!,
            vcat([1.0, 0.0, 0.0, 0.5, 1.0, 0.0], ones(36)),  # placeholder for u0
            [0.0, 1.0],                                      # placeholder for tspan
            (mu,),
        )
        
        # initialize
        new(
            problem,
            method,
            reltol,
            abstol,
        )
    end
end


"""
Propagate initial state `u0` from `tspan[1]` to `tspan[2]`.
The initial state should be given as `u0 = [x,y,z,vx,vy,vz]`.
"""
function propagate(
    propagator::CR3BPPropagator,
    tspan::Tuple{Real,Real},
    u0::Vector;
    callback = nothing,
)
    @assert length(u0) == 6 "u0 should be length 6, not $(length(u0))!"
    @assert length(tspan) == 2 "tspan should be length 2, not $(length(tspan))!"

    # mutate ODEProblem 
    modified_problem = remake(propagator.problem;
                              u0 = u0,
                              tspan = tspan)

    # solve and return results
    return solve(modified_problem,
                 propagator.method;
                 callback = callback,
                 reltol = propagator.reltol,
                 abstol = propagator.abstol)
end


"""
Propagate initial state `u0` from `tspan[1]` to `tspan[2]`.
The initial state should be given as `u0 = [x,y,z,vx,vy,vz]`.
"""
function propagate(
    propagator::CR3BPPropagatorSTM,
    tspan::Tuple{Real,Real},
    u0::Vector;
    callback = nothing,
    stm0 = nothing,
)
    @assert length(u0) == 6 "u0 should be length 6, not $(length(u0))!"
    @assert length(tspan) == 2 "tspan should be length 2, not $(length(tspan))!"

    if isnothing(stm0) == true
        u0_stm = vcat(u0, reshape(I(6), 36))
    else
        u0_stm = vcat(u0, reshape(stm0, 36))
    end

    # mutate ODEProblem 
    modified_problem = remake(propagator.problem;
                              u0 = u0_stm,
                              tspan = tspan)

    # solve and return results
    return solve(modified_problem,
                 propagator.method;
                 callback = callback,
                 reltol = propagator.reltol,
                 abstol = propagator.abstol)
end