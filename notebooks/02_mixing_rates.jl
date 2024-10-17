using DrWatson
using FHist
using HDF5
using JLD2
using DataFrames
using DimensionalData

function load_obs_hist(; file = datadir("obs/") * "wind_hist3d.h5")
    h = h5readhist(file, "hist")
    xs = bincenters(h)
    ps = bincounts(h)
    vs=2 .^ xs[1]
    βs=xs[2] ./ 2
    θs=xs[3]
    return DimArray(ps, (v=vs, β=βs, θ=θs))
end

load_tm_df(; file = datadir() * "/tm_dir=test_alg.jld2") = load(file)["df"]

p = load_obs_hist()
tm_df = load_tm_df()

rotate180(M) = reverse(reverse(M, dims=1), dims=2)

"""
vP is the velocity of the particle, v is the normalized velocity

We find out that for vs smaller than 0.1 and larger than 200, the simulation results are the same as the closest value in the range [0.1, 200]. We filter out these values to save computation time.
"""
function get_tm(v0, β, θ; vP = 512, df = tm_df, vmin = 0.1, vmax = 200)
    v = vP / v0
    if v < vmin
        v = vmin
    elseif v > vmax
        v = vmax
    end

    sdf = @subset(df, :v .== vP/v0 , :β .== β, :θ .== θ)
    if isempty(sdf)
        @info "No data found for v=$v, β=$β, θ=$θ"
        return NaN
    end
    return mean(sdf.tm)
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

tm_stats = mapreduce(+, Iterators.product(vs, βs, θs)) do (v, β, θ)
    p[v=At(v), β=At(β), θ = At(θ)] * get_tm(v, β, θ)
end