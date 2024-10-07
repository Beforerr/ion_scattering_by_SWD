"""
    v_init(v, w, ϕ, r)

Initialize the velocity vector of the particle with magnitude `v`, pitch angle `w`, and azimuthal angle `ϕ` at position `r`.
"""
function v_init(v, w, ϕ, e_para, e_perp1, e_perp2)
    v_para = v * w
    v_perp1 = v * sqrt(1 - w^2) * cos(ϕ)
    v_perp2 = v * sqrt(1 - w^2) * sin(ϕ)
    return v_para * e_para + v_perp1 * e_perp1 + v_perp2 * e_perp2
end

function v_init(v, w, ϕ, B; e1 = ez)
    e_para = B |> normalize  # parallel to B
    e_perp1 = e1 × e_para |> normalize
    e_perp2 = e_para × e_perp1 |> normalize
    return v_init(v, w, ϕ, e_para, e_perp1, e_perp2)
end

v_init(v, w, ϕ, r, B_field; e1 = ez) = v_init(v, w, ϕ, B_field(r); e1 = e1)

# (Optional) Guiding Center Calculation
function solve_gc(param, stateinit, dstate, state, t)
    gc = get_gc(param)
    gc_x0 = gc(stateinit)
    prob_gc = ODEProblem(trace_gc!, gc_x0, tspan, (param..., sol))
    return solve(prob_gc, Tsit5(); save_idxs=[1, 2, 3])
end