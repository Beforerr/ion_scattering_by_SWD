using DrWatson
@quickactivate
using Revise
using CairoMakie
using Beforerr
includet("../src/main.jl")


dir = "simulations"
df = get_result(; dir);

θs = vals(df, :θ)
βs = vals(df, :β)

scale = scales(Color=(;
    lowclip=:transparent,
    colorrange=(0.01, 10),
))

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
let dir = "test"
    df = get_result(; dir);
    ldf, rdf = split_results(df);
    pa_pair_plot(ldf, rdf)
    easy_save("tm/example")
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