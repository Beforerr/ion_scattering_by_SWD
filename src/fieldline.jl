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

function solve_fl(r, B_field; tspan = (0.0, 100.0), alg = Tsit5())
    prob = ODEProblem(field_line_ode!, r, tspan, (B_field,))
    return solve(prob, alg)
end