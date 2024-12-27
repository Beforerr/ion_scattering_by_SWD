using DrWatson
@quickactivate
using Revise
using CairoMakie
using Beforerr
include("../src/main.jl")


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

let figure = (; size=(800, 400))
    pa_layer() = mapping(xyw...) * density_layer()
    subset_outside(df) = @subset(df, outside.(:u1; z_init_0=3.5))

    sdf = @subset(df, :θ .== 85, :β .== 60, :v .!= 8) |> subset_outside
    pa_pair_hist(sdf; figure, layer=pa_layer())
    easy_save("tm/example_subset")
end

let figure = (; size=(400, 800))
    sdf = @subset(df, :θ .== 85, :sign .== 1)
    fig = Figure(; figure)
    layer = data(sdf) * mapping(μ0, ϕ0) * visual(Heatmap) * col_row_mapping(:v)
    grids = sdraw!(fig[1, 1], layer * (:dR_perp_asym,), :v)
    colorbar!(fig[1, end+1], grids[end])
    fig
end

let figure = (; size=(400, 800))
    sdf = @subset(df, :θ .== 85, :sign .== 1, :β .== 90)
    fig = Figure(; figure)
    layer = data(sdf) * mapping(μ0, ϕ0) * visual(Heatmap) * col_row_mapping(:v)
    grids = sdraw!(fig, layer * (:dR_perp_asym_norm,), :v; add_cb=true)
    fig
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