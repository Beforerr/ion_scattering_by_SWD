using StatsBase, Distributions
using DimensionalData
using DataFrames, DataFramesMeta

"""
    EmpiricalDistribution(data::Vector{T} where T <: Real)

Create a discrete empirical distribution based on observations.
"""
function EmpiricalDistribution(data)
    data = sort(data) #sort the observations
    empirical_cdf = ecdf(data) #create empirical cdf
    data_clean = unique(data) #remove duplicates to avoid allunique error
    cdf_data = empirical_cdf.(data_clean) #apply ecdf to data
    pmf_data = vcat(cdf_data[1], diff(cdf_data)) #create pmf from the cdf
    DiscreteNonParametric(data_clean, pmf_data) #define distribution
end

function rand_jump(x0, da)
    d = da[x=Near(x0)]
    rand(d)
end

function rand_jumps!(x::AbstractVector, args...; n=1)
    for _ in 1:n
        x_new = x[end] + rand_jump(x[end], args...)
        push!(x, x_new)
    end
    return x
end

rand_jumps(x0, args...; kw...) = rand_jumps!([x0], args...; kw...)

function ecdf_da(df, x, y)
    ecdfs = combine(
        groupby(df, x),
        y => EmpiricalDistribution;
        renamecols=false
    )
    sort!(ecdfs, x)
    DimArray(ecdfs[!, y], (x=ecdfs[!, x],))
end

function df_rand_jumps(results, x, y; n = 100, counts=1024)
    da = ecdf_da(results, x, y)
    min_x, max_x = extrema(ecdfs[!, x])
    xs = range(min_x, max_x, length=counts)
    stack(rand_jumps.(xs; n, da=da))
end