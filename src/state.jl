function _init_state(v, w, ϕ, B)
    r₀ = init_pos(v)
    v₀ = init_v(v, w, ϕ, B(r₀))
    return [r₀..., v₀...]
end

function init_state(B, v, wϕs::Vector)
    return map(wϕs) do wϕ
        _init_state(v, wϕ..., B)
    end
end

function init_state(B, v, args...; kwargs...)
    wϕs = w_ϕ_pairs(args...; kwargs...)
    return init_state(B, v, wϕs)
end

function local_B_coord(B; e1=ez)
    e_para = B |> normalize  # parallel to B
    e_perp1 = e1 × e_para |> normalize
    e_perp2 = e_para × e_perp1 |> normalize
    return e_para, e_perp1, e_perp2
end

"""
    init_v(v, w, ϕ, ...)

Initialize the velocity vector of the particle with magnitude `v`, pitch angle `w`, and azimuthal angle `ϕ` at position `r`.
"""
function init_v(v, w, ϕ)
    v_para = v * w
    v_perp1 = v * sqrt(1 - w^2) * cos(ϕ)
    v_perp2 = v * sqrt(1 - w^2) * sin(ϕ)
    return [v_para, v_perp1, v_perp2]
end

function init_v(v, w, ϕ, e_para, e_perp1, e_perp2)
    return [e_para e_perp1 e_perp2] * init_v(v, w, ϕ)
end

init_v(v, w, ϕ, B0; e1=ez) = init_v(v, w, ϕ, local_B_coord(B0; e1=e1)...)
init_v(v, w, ϕ, r, B; e1=ez) = init_v(v, w, ϕ, B(r); e1=e1)

"""
Generate a grid of w and phi pairs for the particles.
"""
function w_ϕ_pairs(; Nw=8, Nϕ=8)
    # Generate Grid Points
    ws = range(0, 1, length=Nw) # 0:
    ϕs = range(0, 2π, length=Nϕ + 1)[1:end-1]
    # Flatten the grid to create lists of w0 and phi0 for each particle
    w_all = repeat(ws, inner=Nϕ)
    ϕ_all = repeat(ϕs, outer=Nw)
    return collect(zip(w_all, ϕ_all))
end

"""
Generate a grid of w and phi pairs for the particles.
"""
function w_ϕ_pairs(w; Nϕ=8)
    # Generate Grid Points
    ϕs = range(0, 2π, length=Nϕ + 1)[1:end-1]
    return collect(zip(repeated(w), ϕs))
end

w_ϕ_pairs(w, ϕ) = [(w, ϕ)]