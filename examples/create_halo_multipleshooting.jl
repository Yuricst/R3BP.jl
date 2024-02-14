"""
Example script to create collinear halo with multiple shooting
"""

using GLMakie
using DifferentialEquations

include(joinpath(@__DIR__, "..", "src", "R3BP.jl"))

# create propagator object
params = R3BP.get_cr3bp_param(399, 301)
prop = R3BP.CR3BPPropagator(Tsit5(), params.mu)
prop_stm = R3BP.CR3BPPropagatorSTM(Tsit5(), params.mu)

# create initial guess
Az_km = 12000.0
initial_guess = R3BP.halo_analytical_construct(mu, 2, Az_km, params.lstar, 3)

# solve via mlultiple shooting
n_nodes = 6
indices = round.(Int, range(1, stop=length(initial_guess.ts), length=n_nodes))
tofs = Float64[]
for i = 1:n_nodes - 1
    if i == 1
        push!(tofs, initial_guess.ts[indices[i+1]])
    else
        push!(tofs, initial_guess.ts[indices[i+1]] - initial_guess.ts[indices[i]])
    end
end
x0s = [initial_guess.xs[indices,:][i,:] for i in 1:n_nodes]

lpo_sol = R3BP.multiple_shooting(
    prop_stm,
    x0s,
    tofs,
    fix_time = true,
    periodic = true,
    fix_x0 = false,
    fix_xf = false,
    verbose = true
)
@show lpo_sol.flag

# propagate initial guess & converged solution
sol_ig = R3BP.propagate(prop, (0.0, initial_guess.period), initial_guess.x0)
sol_converged = R3BP.propagate(prop, (0.0, initial_guess.period), lpo_sol.xs[1])

# create Lagrange points for plotting
LPs = R3BP.lagrange_points(mu)

# plot with GLMakie
fig = Figure(size=(600,400), fontsize=22)
ax1 = Axis3(fig[1, 1], aspect=:data)
scatter!(ax1, LPs[2,1], LPs[2,2], LPs[2,3], markersize=10, color=:red)
#lines!(ax1, sol_ig[1,:], sol_ig[2,:], sol_ig[3,:])
lines!(ax1, sol_converged[1,:], sol_converged[2,:], sol_converged[3,:])

display(fig)