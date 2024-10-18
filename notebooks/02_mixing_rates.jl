using DrWatson
using FHist
using HDF5
using JLD2
using DataFrames
using DimensionalData
using StatsBase: mean, midpoints
using Beforerr

BINS = 45

function load_obs_hist(; file=datadir("obs/") * "wind_hist3d.h5")
    h = h5readhist(file, "hist")
    bc = bincenters(h)
    vs = 2 .^ bc[1]
    βs = bc[2] ./ 2
    θs = bc[3]
    return DimArray(bincounts(h), (v=vs, β=βs, θ=θs))
end

w_edges = range(-1, 1, length=BINS + 1)
w_centers = midpoints(w_edges)

tm_file(; bins=BINS, dir="simulations") = datadir() * "/tm_" * savename(@dict bins dir) * ".jld2"
load_tm_df(; kw...) = load(tm_file(; kw...))["df"]

p = load_obs_hist()
tm_df = load_tm_df()

VMIN = 0.25
VMAX = 256

"""
vP is the velocity of the particle, v is the normalized velocity

We find out that for vs smaller than 0.1 and larger than 200, the simulation results are the same as the closest value in the range [0.1, 200]. We filter out these values to save computation time.
"""
function get_tm(v0, β, θ, vP; df=tm_df, vmin=VMIN, vmax=VMAX)
    v = vP / v0
    if v < vmin
        v = vmin
    elseif v > vmax
        v = vmax
    end

    sdf = @subset(df, :v .== v, :β .== β, :θ .== θ)
    tms = sdf[:, :tm]
    if isempty(tms)
        @warn "No data found for v=$v, β=$β, θ=$θ, v0=$v0"
    elseif length(tms) != 1
        @warn "Multiple data found for v=$v, β=$β, θ=$θ, v0=$v0"
    end
    tm = tms[1]
    return mean([tm, rot180(tm)])
end

"""
Check if the parameter range of the transition matrix covers the parameter range of the observation data
"""
function check_coverrange(p, tm_df)
    p_θs = dims(p, :θ)
    p_βs = dims(p, :β)
    p_vs = dims(p, :v)
    sim_θs = tm_df.θ |> unique
    sim_βs = tm_df.β |> unique
    sim_vs = tm_df.v |> unique
    @info p_θs ⊆ sim_θs
    @info p_βs ⊆ sim_βs
    @info p_vs ⊆ sim_vs
end

check_coverrange(p, tm_df)

using Unitful
using Unitful: mp

V_UNIT = u"km/s"

function v2E(v; m =mp)
    E = 0.5m * v^2
    round(Int, u"eV", E)
end

vPs = 4 .^ (3:6)
tm_stats_vPs = map(vPs) do vP
    tm_stats = mapreduce(+, Iterators.product(dims(p)...)) do (v, β, θ)
        p[v=At(v), β=At(β), θ=At(θ)] * get_tm(v, β, θ, vP)
    end
end

tm_stats_vPs_df = mapreduce(vcat, vPs, tm_stats_vPs) do v, tm_stats
    df = DataFrame(DimArray(tm_stats, (w0=w_centers, w1=w_centers)))
    insertcols!(df, :v => v * V_UNIT)
end

using AlgebraOfGraphics

# Plot the transition matrix weighted by the observation data
# As the observation data is dominated by θ->90° (small $B_n$) and β->45°, the allover transition matrix is dominated by these values

vP_map = :v => (v -> "vₚ = $(v) (E = $(v2E(v)))") # particle velocity
plt = data(tm_stats_vPs_df) * mapping(xyw..., :value, layout=vP_map) * visual(Heatmap)
draw(plt; axis=w_axis, colorbar=(; scale=log10))
easy_save("tm_stats_vPs")



tm_stats_das = map([1, 2, 3, 4, 5]) do i
    DimArray(tm_stats^i, (w₀=w_centers, w₁=w_centers))
end
using CairoMakie
f, ax, hm = plot(tm_stats_das[5]; axis=w_axis, colorscale=log10)
Colorbar(f[1, 2], hm)
f