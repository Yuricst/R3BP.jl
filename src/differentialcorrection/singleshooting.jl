"""
Single shooting differential correction algorithm for periodic orbits in R3BP systems
"""



# -------------------------------------------------------------------------------------------- #
# single-shooting differnetial correction
struct Struct_out_ssdc
    x0
    period::Float64
    sol::ODESolution
    flag::Int
    fiters::Vector
end

struct SingleShootingSolution
    x0::Vector
    period::Float64
    flag::Int
    fiters::Vector
end

"""
    ssdc_periodic_xzplane(p, x0, period0; kwargs...)

Single-shooting differential correction for periodic trajectory with symmetry across xz-plane

# Arguments
    p (tuple): parameters for DifferentialEquations
    x0 (Array):
    period0 (float): period of LPO; shooting is done based on perpendicular plane crossing at period0/2
    kwargs:
        maxiter
        reltol
        abstol
        method
        fix
        tolDC
        system (str): "cr3bp" or "er3bp"
        stm_option (str): "analytical" or "ad"

# Returns
    (struct): struct with fields: x0, period, sol, flag, fiters
"""
function single_shooting_xzplane(
    propagator::CR3BPPropagatorSTM,
    x0::Vector,
    period0::Real;
    maxiter::Int = 15,
    fix::String = "period",
    tolDC::Float64 = 1e-12,
    verbose::Bool = false,
)
    # initialize with array and period
    x0iter = deepcopy(x0)
    period = deepcopy(period0)
    dstatef = zeros(length(x0iter))

    # ----- iterate until convergence ----- #
    idx = 1
    flag = 0
    fiters = []
    while idx < maxiter+1
        # propagate state & STM
        sol = R3BP.propagate(propagator, (0.0, period/2), x0iter)
        statef = sol.u[end][1:length(x0iter)]
        stm = reshape(sol.u[end][length(x0iter)+1:end], (length(x0iter),length(x0iter)))';
        propagator.eom!(dstatef, statef, propagator.problem.p, sol.t[end])
        
        # residual vector
        if length(x0iter)==4
            #        y         vx
            ferr = [ statef[2] statef[3] ]'
        elseif length(x0iter)==6
            #        y         vx        vz
            ferr = [ statef[2] statef[4] statef[6] ]'
        end

        # create free parameters vector and DF
        if fix=="vy" && length(x0iter)==4
            #      dy / dx    dvx / d(P/2)
            xi = [ x0iter[1]  period/2 ]'
            df = [ [ stm[2,1] dstatef[2] ];
                   [ stm[3,1] dstatef[3] ] ];

        elseif fix=="z" && length(x0iter)==6
            #      / dx      / dvy     / d(P/2)
            xi = [ x0iter[1] x0iter[5] period/2 ]'
            df = [ [ stm[2,1] stm[2,5] dstatef[2] ];
                   [ stm[4,1] stm[4,5] dstatef[4] ];
                   [ stm[6,1] stm[6,5] dstatef[6] ] ];

        elseif fix=="period" && length(x0iter)==4
            #      dy / dx    dvx / dvy
            xi = [ x0iter[1]  x0iter[4] ]'
            df = [ [ stm[2,1] stm[2,4] ];
                   [ stm[3,1] stm[3,4] ] ];

        elseif fix=="period" && length(x0iter)==6
            #      / dx      / dz      / dvy
            xi = [ x0iter[1] x0iter[3] x0iter[5] ]'
            df = [ [ stm[2,1] stm[2,3] stm[2,5] ];
                   [ stm[4,1] stm[4,3] stm[4,5] ];
                   [ stm[6,1] stm[6,3] stm[6,5] ] ];
        end

        # correct state
        xii = xi - inv(df)*ferr
        push!(fiters, norm(ferr))

        # update state
        if fix=="vy" && length(x0iter)==4
            #      x         P/2
            x0iter[1] = xii[1];
            period    = 2xii[2];
        elseif fix=="z" && length(x0iter)==6
            #      x         vy        P/2
            x0iter[1] = xii[1];
            x0iter[5] = xii[2];
            period    = 2xii[3];
        elseif fix=="period" && length(x0iter)==4
            #      x         vy
            x0iter[1] = xii[1];
            x0iter[4] = xii[2];
        elseif fix=="period" && length(x0iter)==6
            #      x         z         vy
            x0iter[1] = xii[1];
            x0iter[3] = xii[2];
            x0iter[5] = xii[3];
        end

        # check convergence
        if norm(ferr) < tolDC
            if verbose==true
                @printf("Converged at iteration %i\n", idx)
            end
            flag = 1
            break
        else
            if verbose==true
                @printf("Iteration %i: residual: %s\n", idx, norm(ferr))
            end
            idx += 1
        end
    end

    # construct problem
    return SingleShootingSolution(x0iter, period, flag, fiters)
end