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

d = (;θ = 85, β = 47.5, layout=:v)
sdf = @subset(df, :θ .== d.θ, :β .== d.β)
v = histogram(; bins=64)

plt0 = data(sdf) * mapping(xyw...; layout=v_map) * v
fg = draw(plt0; axis=w_axis)
fg

plt = data(subset_leave(sdf)) * mapping(xyw...; layout=v_map) * v
draw(plt; axis=w_axis)
Label(fg.figure[0,1:end], "θ = $(d.θ)°, β = $(d.β)°")
easy_save(savename("pa", d), dir=figure_dir)

plt2 = data(subset_trap(sdf)) * mapping(xyw...; layout=v_map) * v
draw(plt2; axis=w_axis)

layer = data(sdf) * mapping(layout=v_map)
#%%
plt0 = layer * mapping(w0, ϕ0, color=Δα)
fg = draw(plt0)
t = :t1 => "Δt"
plt0 = layer * mapping(w0, ϕ0, color=t)
fg = draw(plt0)

Label(fg.figure[0,1:end], "θ = $(d.θ)°, β = $(d.β)°")