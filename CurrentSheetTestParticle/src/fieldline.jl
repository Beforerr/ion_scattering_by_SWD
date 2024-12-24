"""
Compute the derivative of the position vector along the field line.
drds = B(r) / |B(r)|
"""
function field_line_ode!(drds, r, p, t)
    B_field = p[1]
    B = B_field(r)
    B_norm = norm(B)
    drds[1] = B[1] / B_norm
    drds[2] = B[2] / B_norm
    drds[3] = B[3] / B_norm
end

function solve_fl(r, B_field; tspan=(0.0, 100.0), alg=Tsit5(), kwargs...)
    prob = ODEProblem(field_line_ode!, r, tspan, (B_field,); kwargs...)
    return solve(prob, alg)
end

function field_lines(sol, B; kwargs...)
    gc = get_gc_func(B)
    gc0 = gc(sol[1])
    gcf = gc(sol[end])

    isoutofdomain = (u, p, t) -> abs(u[3]) > maximum(abs.(sol[3, :]))
    tmax = sol.t[end]
    fl0_sol = solve_fl(gc0, B; tspan=(0.0, tmax), isoutofdomain, kwargs...)
    flf_sol = solve_fl(gcf, B; tspan=(0.0, -tmax), isoutofdomain, kwargs...)
    return fl0_sol, flf_sol
end

distance(sol1::ODESolution, sol2::ODESolution) = distance(sol1[1:3, :]', sol2[1:3, :]')

"""field lines distance"""
field_lines_distance(sol, B) = distance(field_lines(sol, B)...)