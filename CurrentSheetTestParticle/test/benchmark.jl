using CurrentSheetTestParticle
using TestParticle
using OrdinaryDiffEq
using DiffEqGPU
using BenchmarkTools

# determine the os and set the backend accordingly

@static if Sys.ismac() 
    using Metal
    backend = Metal.MetalBackend()
end

const DEFAULT_TSPAN = CurrentSheetTestParticle.DEFAULT_TSPAN
const DEFAULT_DIFFEQ_KWARGS = CurrentSheetTestParticle.DEFAULT_DIFFEQ_KWARGS
const DEFAULT_SOLVER = CurrentSheetTestParticle.DEFAULT_SOLVER
const E = CurrentSheetTestParticle.E

v = 8
B = RD_B_field
u0s = init_state(B, v; Nϕ = 64)
isoutofdomain = CurrentSheetTestParticle.isoutofdomain_params(v)

function solve_params_gpu(B, u0s::Vector; E=E, alg=DEFAULT_SOLVER, tspan=DEFAULT_TSPAN, diffeq=DEFAULT_DIFFEQ_KWARGS, kwargs...)
    solve_kwargs = merge(diffeq, kwargs)

    param = prepare(E, B; species=User)
    prob = ODEProblem(trace_normalized!, u0s[1], tspan, param)

    ensemble_prob = EnsembleProblem(prob, u0s; safetycopy=false)
    solve(ensemble_prob, alg, EnsembleGPUArray(); trajectories=length(u0s), solve_kwargs...)
end

save_everystep = false
@btime solve_params_gpu(B, u0s; save_everystep, isoutofdomain)
@btime solve_params(B, u0s; save_everystep, isoutofdomain)