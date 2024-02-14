"""
Example script to create collinear halo
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

# solve via single shooting
lpo_sol = R3BP.single_shooting_xzplane(
    prop_stm,
    initial_guess.x0,
    initial_guess.period,
    fix = "z",
)
@show lpo_sol.flag

# propagate initial guess & converged solution
sol_ig = R3BP.propagate(prop, (0.0, initial_guess.period), initial_guess.x0)
sol_converged = R3BP.propagate(prop, (0.0, lpo_sol.period), lpo_sol.x0)

# create Lagrange points for plotting
LPs = R3BP.lagrange_points(mu)

# plot with GLMakie
fig = Figure(size=(600,400), fontsize=22)
ax1 = Axis3(fig[1, 1], aspect=:data)
scatter!(ax1, LPs[2,1], LPs[2,2], LPs[2,3], markersize=10, color=:red)
#lines!(ax1, sol_ig[1,:], sol_ig[2,:], sol_ig[3,:])
lines!(ax1, sol_converged[1,:], sol_converged[2,:], sol_converged[3,:])

display(fig)