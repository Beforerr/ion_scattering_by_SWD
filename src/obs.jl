using DrWatson
using DimensionalData
using DimensionalData: dims
using HDF5

function load_obs_hist(; file=datadir("obs/") * "wind_hist3d.h5")
    h = h5readhist(file, "hist")
    bc = bincenters(h)
    vs = 2 .^ bc[1]
    βs = bc[2] ./ 2
    θs = bc[3]
    return DimArray(bincounts(h), (v=vs, β=βs, θ=θs))
end

function dargmax(p)
    # Note: iteration is deliberately unsupported for CartesianIndex.
    name.(dims(p)), map(getindex, dims(p), Tuple(argmax(p)))
end

function check_obs_hist(p)
    # check the maxximum value of the observation data
    @info maximum(p)
    @info dargmax(p)
end