using DrWatson
using DataFrames
using DataFramesMeta
using CurrentSheetTestParticle

include("io.jl")
include("plot.jl")

subset_v(df) = @subset(df, :v .<= 128)
subset_v(v::Number) = df -> @subset(df, :v .<= v)
subset_β(itr) = df -> @rsubset(df, :β ∈ itr)

#%%
vals(df, s) = unique(df[!, s]) |> sort