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

d = (;θ = 85, β = 47.5, layout = :v)
sdf = @subset(df, :θ .== d.θ, :β .== d.β)
v = histogram(; bins=64)
plt = data(subset_leave(sdf)) * mapping(xyw...; layout=v_map) * v
plt2 = data(subset_trap(sdf)) * mapping(xyw...; layout=v_map) * v
draw(plt; axis=w_axis)
easy_save(savename("pa", d), dir=figure_dir)