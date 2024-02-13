"""
Example script to propagate state in CR3BP
"""

using GLMakie
using DifferentialEquations

include(joinpath(@__DIR__, "..", "src", "R3BP.jl"))

# create propagator object
mu = 0.01215058426994
prop = R3BP.CR3BPPropagator(Tsit5(), mu)

# propagate 
u0 = [1.176924090973164, 0.0, -0.060210863312217, 0.0, -0.173836346247689, 0.0];
tspan = (0.0, 3.385326412831325)
sol = R3BP.propagate(prop, tspan, u0)

# create Lagrange points for plotting
LPs = R3BP.lagrange_points(mu)

# plot with GLMakie
fig = Figure(size=(600,400), fontsize=22)
ax1 = Axis3(fig[1, 1], aspect=:data)
scatter!(ax1, LPs[1:2,1], LPs[1:2,2], LPs[1:2,3], markersize=10, color=:red)
lines!(ax1, sol[1,:], sol[2,:], sol[3,:])

# Generate points on the sphere
nsph = 10
θ = range(0, stop=2π, length=nsph)
ϕ = range(0, stop=π, length=nsph)
R = 1737.4/384400      # approximate size of the Moon, in LU
center = [1-mu, 0, 0]
xsphere = [center[1] + R * cos(θ[i]) * sin(ϕ[j]) for j in 1:nsph, i in 1:nsph]
ysphere = [center[2] + R * sin(θ[i]) * sin(ϕ[j]) for j in 1:nsph, i in 1:nsph]
zsphere = [center[3] + R * cos(ϕ[j]) for j in 1:nsph, i in 1:nsph]
wireframe!(ax1, xsphere, ysphere, zsphere, color=:black, linewidth=1.0)

display(fig)