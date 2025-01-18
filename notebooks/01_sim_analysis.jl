begin
    using DrWatson
    @quickactivate
    using Revise
    using Beforerr
    using Statistics
    include("../src/main.jl")
end

struct vθβ
    v
    θ
    β
end

rename_func(t::vθβ) = L"v_p = %$(t.v) v_0,\ θ = %$(t.θ)^∘,\ β = %$(t.β)^∘"

# ----
# Simulations
# ----
dir = "simulations"
df = get_result(; dir);

θs = vals(df, :θ)
βs = vals(df, :β)

d = (; θ=85, β=47.5, layout=:v)
sdf = @subset(df, :θ .== d.θ, :β .== d.β)
let v = histogram(; bins=64, normalization=:pdf)
    plt0 = data(sdf) * mapping(xyw...; layout=v_map) * v
    draw(plt0, scale; axis=w_axis)
    easy_save(savename("pa", d))
end


# ----
# Example
# ----
dir = "example"
df = get_result(; dir)

let figure = (; size=(72 * 6.5, 72 * 6), figure_padding=0), axis = w_axis, colorbar = (; scale=colorscale)
    layer = mapping(xyw...) * density_layer()
    subset_outside(df) = @subset(df, outside.(:u1; z_init_0=3))
    sdf = @subset(df, +(:v .== 8, :β .== 50, :θ .== 85) .>= 2) |> subset_outside
    sdf.group = vθβ.(sdf.v, sdf.θ, sdf.β)
    scales = tm_scale(; Layout=(; palette=[(2, 2), (2, 1), (1, 1), (1, 2)]))

    draw(data(sdf) * layer * mapping(layout=:group => rename_func), scales; figure, axis, colorbar)
    easy_save("tm/example_subset")
end

# ----
# Test
# ----
colorscale = log10

dir = "test"
df = get_result(; dir)

let
    ldf, rdf = split_results(df)
    pa_pair_hist(ldf, rdf)
    easy_save("tm/example")
end

let figure = (; size=(400, 800)),
    sdf = @subset(df, :θ .== 85, :sign .== 1)

    fig = Figure(; figure)
    layer = data(sdf) * mapping(μ0, ϕ0) * visual(Heatmap) * col_row_mapping(:v)
    grids = sdraw!(fig[1, 1], layer * (:dR_perp_asym,), :v)
    colorbar!(fig[1, end+1], grids[end])
    fig
end

let figure = (; size=(1200, 400)), colorrange = (0, 40)
    sdf = @subset(df, :v .> 1, :θ .> 60, :β .>= 60)
    fig = Figure(; figure)
    layer = data(sdf) * mapping(μ0, ϕ0) * visual(Heatmap) * col_row_mapping(:v)
    grids = sdraw!(fig, layer * (dR_perp_asym_norm,), v_map; scales=scales(Color=(; colorrange)))
    colorbar!(fig[1, end+1], grids[1])
    easy_save("diffusion/dR_perp_asym_norm")
end

# Group by (v, θ, β) and get the average dR_perp_asym_norm
let cols = [:v, :θ, :β]
    sdf = @subset(df, :v .> 1)
    tdf = combine(groupby(sdf, cols), :dR_perp_asym_norm => mean; renamecols=false)
    draw(data(tdf) * mapping(:β, dR_perp_asym_norm; color=v_map, col=θ_map))
    easy_save("diffusion/dR_perp_asym_norm_avg")
end
# ----
# Test different field configurations
# ----
let dir = "test_TD", figure = (; size=(500, 1000))
    df = get_result(; dir)
    pa_layer() = mapping(xyw..., col=:Bfn => string, row=v_map) * density_layer()
    draw(data(df) * pa_layer(), tm_scale(); figure)
    easy_save("tm/example_td")
end

plt = data(subset_leave(sdf)) * mapping(xyw...; layout=v_map) * v
draw(plt; axis=w_axis)
Label(fg.figure[0, 1:end], "θ = $(d.θ)°, β = $(d.β)°")


plt2 = data(subset_trap(sdf)) * mapping(xyw...; layout=v_map) * v
draw(plt2; axis=w_axis)

layer = data(sdf) * mapping(layout=v_map)
#%%
plt0 = layer * mapping(w0, ϕ0, color=Δα)
fg = draw(plt0)
t = :t1 => "Δt"
plt0 = layer * mapping(w0, ϕ0, color=t)
fg = draw(plt0)

Label(fg.figure[0, 1:end], "θ = $(d.θ)°, β = $(d.β)°")