using DataFramesMeta
using CurrentSheetTestParticle
using CurrentSheetTestParticle: B, DEFAULT_Z_INIT_0, DEFAULT_TSPAN, DEFAULT_SIGN
include("utils.jl")

#%%
outside(u; z_init_0=DEFAULT_Z_INIT_0) = abs(u[3]) > z_init_0
subset_leave(df) = @subset(df, :t1 .!= :tmax)
subset_outside(df) = @subset(df, outside.(:u1))
subset_trap(df) = @subset(df, :t1 .== :tmax)

function get_result_dfs(; dir="simulations")
    df = collect_results(datadir(dir))
    cols = names(df)
    "tspan" ∈ cols ? @transform!(df, :tmax = last.(:tspan)) : insertcols!(df, :tmax => DEFAULT_TSPAN[2])
    "sign" ∈ cols || insertcols!(df, :sign => DEFAULT_SIGN)
    "Bfn" ∈ cols || insertcols!(df, :Bfn => RotationDiscontinuity)
    @chain df begin
        @transform!(:sign = coalesce.(:sign, DEFAULT_SIGN))
        @rtransform!(:B = B(:Bfn(θ=:θ, β=:β, sign=:sign)))
    end
end

function get_result(; dir="simulations")
    dfs = get_result_dfs(; dir)
    results = map(get_result, eachrow(dfs))
    vcat(results...)
end

function get_result!(r::DataFrameRow, params=[:θ, :β, :v, :tmax, :B, :alg, :sign, :Bfn])
    params = string.(params) ∩ names(r)
    @chain r[:results] begin
        insertcols!(Pair.(params, Array(r[params]))...; makeunique=true)
        process_result!
    end
end

function split_results(df)
    ldf = @subset(df, :sign .== 1)
    rdf = @subset(df, :sign .== -1)
    return ldf, rdf
end