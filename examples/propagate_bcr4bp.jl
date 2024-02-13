"""
Example script to propagate state in CR3BP
"""

using GLMakie
using DifferentialEquations

include(joinpath(@__DIR__, "..", "src", "R3BP.jl"))

# create propagator object
params = R3BP.get_bcr4bp_param("399", "301")
prop = R3BP.BCR4BPPropagator(Tsit5(), params.mu, params.μ_3, params.a, params.ω_s)

# propagate 
t0 = 0.0
u0 = [0.9853006895390061, 0.0021723003109550003, 0.0034046027917,
      -1.198746524193054, 1.038423839870351, -1.559959799100589]
tspan = (0.0, -22.970560065381424)
sol = R3BP.propagate(prop, t0, tspan, u0)
@show sol.u[end]

# create Lagrange points for plotting
LPs = R3BP.lagrange_points(mu)

# plot with GLMakie
fig = Figure(size=(600,400), fontsize=22)
ax1 = Axis3(fig[1, 1], aspect=:data)
scatter!(ax1, LPs[1:5,1], LPs[1:5,2], LPs[1:5,3], markersize=10, color=:red)
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