using TestParticle
using TestParticle: c, qᵢ, mᵢ
using LinearAlgebra
using OrdinaryDiffEq
using StaticArrays
using CairoMakie
using DrWatson
using SciMLBase: AbstractSciMLSolution

include("../src/main.jl")

m = mᵢ

B₀ = 1e-8             # Magnetic field [T]
l₀ = 6e8               # Length [m]
Ω = abs(qᵢ) * B₀ / mᵢ # [1/s]
t₀ = 1 / Ω

# Step 1: Define the Magnetic Field B(x)
z_init = -16

# Simulation Parameters
alg = Vern9()
tspan = (0.0, 512)
Nw = 10           # Number of w points
Nϕ = 8           # Number of phi points
# Generate Grid Points
ws = range(0, 1, length=Nw) # 0:
ϕs = range(0, 2π, length=Nϕ + 1)[1:end-1]
# Flatten the grid to create lists of w0 and phi0 for each particle
w0_all = repeat(ws, inner=Nϕ)
phi0_all = repeat(ϕs, outer=Nw)
# Total number of particles
num_particles = Nw * Nϕ

function RD_B_field(r, α, β, B=1)
    z = r[3]
    ψ = β * tanh(z)
    Bz = B * cos(α)
    Bx = B * sin(α) * cos(ψ)
    By = B * sin(α) * sin(ψ)
    return SVector(Bx, By, Bz)
end

# Step 3: Define the Electric Field (if any)
E_field(x) = SVector(0.0, 0.0, 0.0)

ez = [0, 0, 1]

function isoutofdomain(u, p, t)
    z = u[3]
    return abs(z) > 1.5 * abs(z_init) ? true : false
end

"Set initial state for EnsembleProblem."
function prob_func(prob, i, repeat)
    # Compute the initial velocity for the i-th ensemble member
    v0 = v_init(v, w0_all[i], phi0_all[i], r₀)
    prob = @views remake(prob, u0=[prob.u0[1:3]..., v0...])
end

function sim(α, β, v; solve_field=false)
    B_field = r -> RD_B_field(r, α, β)

    # Initial Phase Space (Position is common; Velocity will be set per ensemble member)
    r₀ = [0, 0, z_init]

    vs = map(w0_all, phi0_all) do w, ϕ
        v_init(v, w, ϕ, r₀, B_field)
    end

    u0s = map(vs) do v
        [r₀..., v...]
    end

    # Prepare the simulation parameters
    param = prepare(E_field, B_field; species=User)
    prob = ODEProblem(trace_normalized!, u0s[1], tspan, param)
    ensemble_prob = EnsembleProblem(prob, u0s; safetycopy=false)
    # Solve the ODE
    sols = solve(ensemble_prob, alg, EnsembleThreads(); trajectories=num_particles, isoutofdomain)
    if !solve_field
        return sols
    else
        sol_field = solve_fl(r₀, B_field)
        return sols, sol_field
    end
end

function makesim(d::Dict)
    @unpack α, β, v = d
    fulld = copy(d)
    result = sim(α, β, v)
    return Dict(fulld..., :result => result)
end

dicts = [
    Dict(:α => π / 4, :β => π / 2, :v => 6),
    Dict(:α => π / 2, :β => π / 2, :v => 6),
    Dict(:α => π / 2 - 0.1, :β => π / 2, :v => 0.1),
]


for d in dicts
    f = makesim(d)
    wsave(datadir("simulations", savename(d, "jld2")), f)
end