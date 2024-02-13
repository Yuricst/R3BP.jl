"""
Propagator classes with BCR4BP
"""


"""
Propagator for BCR4BP
"""
mutable struct BCR4BPPropagator <: R3BPPropagator
    problem
    parameters
    method
    reltol::Float64
    abstol::Float64

    function BCR4BPPropagator(
        method,
        mu::Real,
        mu3::Real,
        as::Real,
        ωs::Real;
        reltol::Float64=1e-12,
        abstol::Float64=1e-12
    )
        # construct ODE problem
        parameters = [mu, mu3, 0.0, as, ωs]
        problem = ODEProblem(
            rhs_bcr4bp_sv!,
            [1.0, 0.0, 0.0, 0.5, 1.0, 0.0],  # placeholder for u0
            [0.0, 1.0],                      # placeholder for tspan
            parameters,                      # placeholder used for t0
        )

        # initialize
        new(
            problem,
            parameters,
            method,
            reltol,
            abstol,
        )
    end
end


"""
Propagator for BCR4BP with STM
"""
mutable struct BCR4BPPropagatorSTM <: R3BPPropagator
    problem
    parameters
    method
    reltol::Float64
    abstol::Float64

    function BCR4BPPropagator(
        method,
        mu::Real,
        mu3::Real,
        as::Real,
        ωs::Real;
        reltol::Float64=1e-12,
        abstol::Float64=1e-12
    )
        # construct ODE problem
        parameters = [mu, mu3, 0.0, as, ωs]
        problem = ODEProblem(
            rhs_bcr4bp_svstm!,
            vcat([1.0, 0.0, 0.0, 0.5, 1.0, 0.0], ones(36)),  # placeholder for u0
            [0.0, 1.0],                                      # placeholder for tspan
            parameters,                                      # placeholder used for t0
        )
        
        # initialize
        new(
            problem,
            parameters,
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
    propagator::BCR4BPPropagator,
    t0::Real,
    tspan::Tuple{Real,Real},
    u0::Vector;
    callback = nothing,
)
    @assert length(u0) == 6 "u0 should be length 6, not $(length(u0))!"
    @assert length(tspan) == 2 "tspan should be length 2, not $(length(tspan))!"

    # mutate ODEProblem
    propagator.parameters[3] = t0
    modified_problem = remake(propagator.problem;
                              u0 = u0,
                              tspan = tspan,
                              p = propagator.parameters)

    # solve and return results
    @show modified_problem.p
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
    propagator::BCR4BPPropagatorSTM,
    t0::Real,
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
    problem.parameters[3] = t0
    modified_problem = remake(propagator.problem;
                              u0 = u0_stm,
                              tspan = tspan,
                              p = problem.parameters)

    # solve and return results
    return solve(modified_problem,
                 propagator.method;
                 callback = callback,
                 reltol = propagator.reltol,
                 abstol = propagator.abstol)
end