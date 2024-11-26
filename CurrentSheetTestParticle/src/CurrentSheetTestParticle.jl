module CurrentSheetTestParticle
using TestParticle
using OrdinaryDiffEq
using LinearAlgebra
using StaticArrays
using UnPack
using Moshi.Match: @match

export RD_B_field, TD_B_field
export solve_params, dsolve_params
export w_ϕ_pairs, init_state, init_states, init_states_pm, filter_wϕs!
export isoutofdomain_params
export ProblemParams
export trace_normalized_B!, trace_normalized_B

include("field.jl")
include("state.jl")
include("fieldline.jl")
include("pa.jl")
include("equations.jl")

abstol = 1e-7 # Defaults to 1e-5
reltol = 1e-7 # Defaults to 1e-3
maxiters = 1e6 # Defaults to 1e5
dtmin = 1e-4 # Reduce the computation for domain checking (as `isoutofdomain` will reject and reducd time step until a step is accepted)
const DEFAULT_SOLVER = AutoVern9(Rodas4P())
const DEFAULT_DIFFEQ_KWARGS = (; abstol, reltol, maxiters, dtmin)
const DEFAULT_BORIS_KWARGS = (; dt=1e-2, savestepinterval=1)
const DEFAULT_TSPAN = (0, 256)
const ez = SA[0, 0, 1]

@kwdef struct ProblemParams
    Bfn = RD_B_field
    θ = DEFAULT_θ
    β = DEFAULT_β
    sign = DEFAULT_SIGN
    v = 1
    alg = :AutoVern9
    init_kwargs = (; Nw=8, Nϕ=8)
    tspan = DEFAULT_TSPAN
    diffeq = DEFAULT_DIFFEQ_KWARGS
end

@kwdef struct BParams
    Bfn = RD_B_field
    θ = DEFAULT_θ
    β = DEFAULT_β
    sign = DEFAULT_SIGN
end

# Step 3: Define the Electric Field (if any)
const E0 = SVector(0.0, 0.0, 0.0)
E(x) = E0

isoutofdomain_z(z_max) = (u, p, t) -> abs(u[3]) > abs(z_max)

function isoutofdomain_params(v)
    z_init = init_z_pos(v)
    z_max = 2 * z_init
    return isoutofdomain_z(z_max)
end

function _alg(alg)
    @match alg begin
        :AutoVern9 => AutoVern9(Rodas4P())
        :Boris => @warn "Boris algorithm is not implemented yet"
        _ => alg
    end
end

isinplace(f) = first(methods(f)).nargs == 5

"""
Solve the system of ODEs.
"""
function solve_params(B, u0s::AbstractVector; f=trace_normalized_B, E=E, alg=DEFAULT_SOLVER, tspan=DEFAULT_TSPAN, diffeq=DEFAULT_DIFFEQ_KWARGS, kwargs...)
    solve_kwargs = merge(diffeq, kwargs)
    alg = _alg(alg)

    u0s = isinplace(f) ? u0s : [SVector{6}(u0) for u0 in u0s]

    param = prepare(E, B; species=User)
    prob = ODEProblem(f, u0s[1], tspan, param)

    prob_func = (prob, i, repeat=nothing) -> remake(prob, u0=u0s[i])
    ensemble_prob = EnsembleProblem(prob; prob_func, safetycopy=false)
    solve(ensemble_prob, alg, EnsembleThreads(); trajectories=length(u0s), solve_kwargs...)
end

function solve_params_boris(B, u0s::AbstractVector; E=E, tspan=DEFAULT_TSPAN, diffeq=DEFAULT_BORIS_KWARGS, kwargs...)
    solve_kwargs = merge(diffeq, kwargs)

    param = prepare(E, B; species=User)
    prob_func = (prob, i, repeat=nothing) -> remake(prob, u0=u0s[i])
    prob = TraceProblem(u0s[1], tspan, param; prob_func)

    TestParticle.solve(prob, EnsembleThreads(); trajectories=length(u0s), solve_kwargs...)
end

function solve_params(d; kwargs...)
    @unpack θ, β, v, sign, alg, init_kwargs, diffeq, tspan, Bfn = d
    B = Bfn(; θ, β, sign)
    u0s, wϕs = init_states_pm(B, v; init_kwargs...)

    isoutofdomain = isoutofdomain_params(v)
    sol = solve_params(B, u0s; alg, tspan, diffeq, isoutofdomain, kwargs...)
    return sol, (wϕs, B)
end

function solve_params(B, v::Number, args...; init_kwargs=(;), kwargs...)
    u0s, wϕs = init_states_pm(B, v, args...; init_kwargs...)
    isoutofdomain = isoutofdomain_params(v)
    sol = solve_params(B, u0s; isoutofdomain, kwargs...)
    return sol, (wϕs,)
end
end
