using CurrentSheetTestParticle
using BenchmarkTools

save_everystep = false
verbose = false
dtmax = 1e-1
diffeq = (; save_everystep, verbose, dtmax)

p = ProblemParams(
    θ=85, β=45, v=64.0, init_kwargs=(; Nw=32, Nϕ=2)
)

@benchmark solve_params($p; f = trace_normalized_B!, diffeq...) seconds=1
@benchmark solve_params($p; f = trace_normalized_B, diffeq...) seconds=1