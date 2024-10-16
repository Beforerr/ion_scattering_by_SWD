module CurrentSheetTestParticle
using TestParticle
using OrdinaryDiffEq
using LinearAlgebra
using StaticArrays

export RD_B_field
export solve_params, dsolve_params
export w_ϕ_pairs, init_state, init_states, init_states_pm, filter_wϕs!
export isoutofdomain_params

include("field.jl")
include("state.jl")
include("fieldline.jl")
include("pa.jl")

abstol = 1e-7 # Defaults to 1e-5
reltol = 1e-7 # Defaults to 1e-3
maxiters = 1e6 # Defaults to 1e5
const DEFAULT_SOLVER = AutoVern9(Rodas4P())
const DEFAULT_DIFFEQ_KWARGS = (; abstol, reltol, maxiters)
const DEFAULT_BORIS_KWARGS = (; dt=1e-2, savestepinterval=1)
const DEFAULT_TSPAN = (0, 256)
ez = [0, 0, 1]

# Step 3: Define the Electric Field (if any)
const E0 = SVector(0.0, 0.0, 0.0)
E(x) = E0

isoutofdomain_z(u, p, t, z_max) = abs(u[3]) > abs(z_max) ? true : false
isoutofdomain_z(z_max) = (u, p, t) -> isoutofdomain_z(u, p, t, z_max)

function isoutofdomain_params(v)
    z_init = init_z_pos(v)
    z_max = 2 * z_init
    return isoutofdomain_z(z_max)
end

"""
Solve the system of ODEs.

Notes: `v` is needed here for domain checking.
"""
function solve_params(B, u0s::Vector; E=E, alg=DEFAULT_SOLVER, tspan=DEFAULT_TSPAN, diffeq=DEFAULT_DIFFEQ_KWARGS, kwargs...)
    solve_kwargs = merge(diffeq, kwargs)

    param = prepare(E, B; species=User)
    prob = ODEProblem(trace_normalized!, u0s[1], tspan, param)

    ensemble_prob = EnsembleProblem(prob, u0s)
    solve(ensemble_prob, alg, EnsembleThreads(); trajectories=length(u0s), solve_kwargs...)
end

function solve_params(B, v, args...; init_kwargs=(;), kwargs...)
    u0s = init_state(B, v, args...; init_kwargs...)
    solve_params(B, u0s; kwargs...)
end

function solve_params_boris(B, u0s::Vector; E=E, tspan=DEFAULT_TSPAN, diffeq=DEFAULT_BORIS_KWARGS, kwargs...)
    solve_kwargs = merge(diffeq, kwargs)

    param = prepare(E, B; species=User)
    prob_func = (prob, i, repeat = nothing) -> remake(prob, u0 = u0s[i])
    prob = TraceProblem(u0s[1], tspan, param; prob_func)

    TestParticle.solve(prob, EnsembleThreads(); trajectories=length(u0s), solve_kwargs...)
end
end
