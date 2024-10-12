using TestParticle
using DynamicalSystems
using OrdinaryDiffEq
using LinearAlgebra
using StaticArrays
export RD_B_field
export solve_params, dsolve_params
export w_ϕ_pairs, init_state

include("state.jl")
include("fieldline.jl")
include("plot.jl")
include("utils.jl")
include("pa.jl")

alg = Vern9()
abstol = 1e-7 # Defaults to 1e-5
reltol = 1e-7 # Defaults to 1e-3
diffeq = (; alg, abstol, reltol)
tspan = (0, 256)
ez = [0, 0, 1]

"""
Rotating magnetic field. 

# Notes
Left-handed rotating with `sign=1`.
Right-handed rotating with `sign=-1`
"""
function RD_B_field(r, α, β; B=1, sign=1)
    z = r[3]
    ψ = β * tanh(z)
    Bz = B * cos(α)
    Bx = B * sin(α) * sin(ψ)
    By = sign * B * sin(α) * cos(ψ)
    return SVector(Bx, By, Bz)
end

RD_B_field(; B=1, α, β, sign) = r -> RD_B_field(r, α, β; B, sign)

# Step 3: Define the Electric Field (if any)
E_field(x) = SVector(0.0, 0.0, 0.0)

# TODO: check the initial position effect
init_z_pos(v, z_init_0=16) = -abs(z_init_0) - 2 * abs(v)
init_pos(v, kwargs...) = [0, 0, init_z_pos(v, kwargs...)]

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
function solve_params(B, v, u0s::Vector; alg=alg, reltol=reltol, kwargs...)
    param = prepare(E_field, B; species=User)
    prob = ODEProblem(trace_normalized!, u0s[1], tspan, param)

    isoutofdomain = isoutofdomain_params(v)

    ensemble_prob = EnsembleProblem(prob, u0s; safetycopy=false)
    solve(ensemble_prob, alg, EnsembleThreads(); trajectories=length(u0s), isoutofdomain, reltol, kwargs...)
end

function solve_params(B, v, args...; init_kwargs=(;), kwargs...)
    u0s = init_state(B, v, args...; init_kwargs...)
    solve_params(B, v, u0s; kwargs...)
end

"""
Using DynamicalSystems.jl to solve the system of ODEs.
"""
function dsolve_params(B, v, args...; alg=alg, saveat=0.25, reltol=reltol)
    u0s = init_state(B, v, args...; init_kwargs...)
    param = prepare(E_field, B; species=User)

    isoutofdomain = isoutofdomain_params(v)
    diffeq = (; alg, isoutofdomain, reltol)

    @Threads.threads for u0 in u0s
        ds = CoupledODEs(trace_normalized!, u0, param; diffeq)
        step!(ds, tspan[2])
    end
end