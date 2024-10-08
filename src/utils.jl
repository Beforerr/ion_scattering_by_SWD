function local_B_coord(B; e1=ez)
    e_para = B |> normalize  # parallel to B
    e_perp1 = e1 × e_para |> normalize
    e_perp2 = e_para × e_perp1 |> normalize
    return e_para, e_perp1, e_perp2
end

"""
    v_init(v, w, ϕ, ...)

Initialize the velocity vector of the particle with magnitude `v`, pitch angle `w`, and azimuthal angle `ϕ` at position `r`.
"""
function v_init(v, w, ϕ)
    v_para = v * w
    v_perp1 = v * sqrt(1 - w^2) * cos(ϕ)
    v_perp2 = v * sqrt(1 - w^2) * sin(ϕ)
    return [v_para, v_perp1, v_perp2]
end

function v_init(v, w, ϕ, e_para, e_perp1, e_perp2)
    return [e_para e_perp1 e_perp2] * v_init(v, w, ϕ)
end

v_init(v, w, ϕ, B; e1=ez) = v_init(v, w, ϕ, local_B_coord(B; e1=e1)...)
v_init(v, w, ϕ, r, B_field; e1=ez) = v_init(v, w, ϕ, B_field(r); e1=e1)

"""
Generate a grid of w and phi pairs for the particles.
"""
function w_ϕ_pairs(Nw=100, Nϕ=100)
    # Generate Grid Points
    ws = range(0, 1, length=Nw) # 0:
    ϕs = range(0, 2π, length=Nϕ + 1)[1:end-1]
    # Flatten the grid to create lists of w0 and phi0 for each particle
    w_all = repeat(ws, inner=Nϕ)
    ϕ_all = repeat(ϕs, outer=Nw)
    return collect(zip(w_all, ϕ_all))
end


# (Optional) Guiding Center Calculation
function solve_gc(param, stateinit, dstate, state, t)
    gc = get_gc(param)
    gc_x0 = gc(stateinit)
    prob_gc = ODEProblem(trace_gc!, gc_x0, tspan, (param..., sol))
    return solve(prob_gc, Tsit5(); save_idxs=[1, 2, 3])
end