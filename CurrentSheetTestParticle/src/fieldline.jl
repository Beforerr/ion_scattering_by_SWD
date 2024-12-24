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

function solve_fl(r, B_field, args...; tspan=DEFAULT_TSPAN, kwargs...)
    prob = ODEProblem(field_line_ode!, r, tspan, (B_field,); kwargs...)
    return solve(prob, args...)
end

function B_field_line(Bf::Function)
    return r -> (B = Bf(r); sign(B[3]) * B)
end

function field_lines(sol, B; kwargs...)
    gc = get_gc_func(B)
    gc0 = gc(sol[1])
    gcf = gc(sol[end])

    isoutofdomain = (u, p, t) -> abs(u[3]) > maximum(abs.(sol[3, :]))
    tmax = 100 * sol.t[end]
    fl0_sol = solve_fl(gc0, B_field_line(B); tspan=(0.0, tmax), isoutofdomain, kwargs...)
    flf_sol = solve_fl(gcf, B_field_line(B); tspan=(0.0, -tmax), isoutofdomain, kwargs...)
    return fl0_sol, flf_sol
end

"""
Distance of two line solutions
"""
distance(sol1::ODESolution, sol2::ODESolution) = distance(sol1[1:3, :]', sol2[1:3, :]')

"""
Asymptotic distance between two line solutions
"""
function field_lines_asym_distance(sol1, sol2, B)
    direction = B([0, 0, Inf])
    p1 = sol1[argmax(sol1[3, :])]
    p2 = sol2[argmax(sol2[3, :])]
    distance(p1, p2, direction)
end

field_lines_asym_distance(sol, B) =
    field_lines_asym_distance(field_lines(sol, B)..., B)

"""
Field lines distances
"""
function field_lines_distance(sol, B)
    fl0_sol, flf_sol = field_lines(sol, B)
    dR_perp_min = distance(fl0_sol, flf_sol)
    dR_perp_asym = field_lines_asym_distance(sol, B)
    return (; dR_perp_min, dR_perp_asym)
end