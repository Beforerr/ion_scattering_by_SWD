using Base.Iterators

Eₖ(u) = 0.5 * norm(u[4:6])^2

# (Optional) Guiding Center Calculation
function solve_gc(param, stateinit, dstate, state, t)
    gc = get_gc(param)
    gc_x0 = gc(stateinit)
    prob_gc = ODEProblem(trace_gc!, gc_x0, tspan, (param..., sol))
    return solve(prob_gc, Tsit5(); save_idxs=[1, 2, 3])
end