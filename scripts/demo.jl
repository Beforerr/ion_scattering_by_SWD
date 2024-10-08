using TestParticle
using TestParticle: c, qᵢ, mᵢ
using LinearAlgebra
using OrdinaryDiffEq
using StaticArrays
using DrWatson
using DataFrames
using Logging

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

# Step 3: Define the Electric Field (if any)
E_field(x) = SVector(0.0, 0.0, 0.0)

ez = [0, 0, 1]

function isoutofdomain(u, p, t)
    z = u[3]
    return abs(z) > 1.5 * abs(z_init) ? true : false
end

function sim(
    α, β, v;
    Nw=100, Nϕ=100,
    solve_field=false
)
    B_field = r -> RD_B_field(r, α, β)

    # Initial Phase Space (Position is common; Velocity will be set per ensemble member)
    r₀ = [0, 0, z_init]

    wϕs = w_ϕ_pairs(Nw, Nϕ)

    vs = map(wϕs) do wϕ
        v_init(v, wϕ..., r₀, B_field)
    end

    u0s = map(vs) do v
        [r₀..., v...]
    end

    # Prepare the simulation parameters
    param = prepare(E_field, B_field; species=User)
    prob = ODEProblem(trace_normalized!, u0s[1], tspan, param)
    ensemble_prob = EnsembleProblem(prob, u0s; safetycopy=false)
    # Solve the ODE
    sols = solve(ensemble_prob, alg, EnsembleThreads(); trajectories=length(u0s), isoutofdomain)
    if !solve_field
        return sols
    else
        sol_field = solve_fl(r₀, B_field)
        return sols, sol_field
    end
end

function sim(d::Dict)
    @unpack α, β, v = d
    return sim(α, β, v)
end


svec(x) = SVector(x...)

function extract_info(sol; add_info=(;))
    u0 = sol.prob.u0 |> svec
    t = sol.t |> svec
    u = sol.u
    u = map(svec, u)

    return Dict(:u0 => u0, :u => u, :t => t, add_info...)
end

function makesim(d::Dict)
    sol = sim(d)
    results = extract_info.(sol.u; add_info=d)
    df = DataFrame(results)
    return Dict(d..., :result => df)
end


Logging.disable_logging(Logging.Warn)

dicts = [
    Dict(:α => π / 4, :β => π / 2, :v => 1),
    Dict(:α => π / 4, :β => π / 2, :v => 8),
    Dict(:α => π / 4, :β => π / 2, :v => 64),
    Dict(:α => π / 2 - 0.1, :β => π / 2, :v => 1),
    Dict(:α => π / 2 - 0.1, :β => π / 2, :v => 8),
    Dict(:α => π / 2 - 0.1, :β => π / 2, :v => 64),
]


for d in dicts
    path = datadir("simulations")
    produce_or_load(makesim, d, path; loadfile=false)
end