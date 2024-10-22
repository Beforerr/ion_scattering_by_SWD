using DataFramesMeta
using CurrentSheetTestParticle
using CurrentSheetTestParticle: DEFAULT_Z_INIT_0, DEFAULT_TSPAN, DEFAULT_SIGN

#%%
outside(u; z_init_0=DEFAULT_Z_INIT_0) = abs(u[3]) > z_init_0
subset_leave(df) = @subset(df, :t1 .!= :tmax)
subset_outside(df) = @subset(df, outside.(:u1))
subset_trap(df) = @subset(df, :t1 .== :tmax)

function get_result_dfs(; dir="simulations")
    df = collect_results(datadir(dir))

    "tspan" ∈ names(df) ? @transform!(df, :tmax = last.(:tspan)) : insertcols!(df, :tmax => DEFAULT_TSPAN[2])
    "sign" ∈ names(df) || insertcols!(df, :sign => DEFAULT_SIGN)

    @rtransform!(df, :B = RD_B_field(; θ=:θ, β=:β, sign=:sign))
end

function get_result(; dir="simulations")
    dfs = get_result_dfs(; dir)
    results = map(get_result, eachrow(dfs))
    vcat(results...)
end

function process_result!(df::AbstractDataFrame)
    @chain df begin
        @rtransform!(:w0 = :wϕ0[1], :ϕ0 = :wϕ0[2], :w1 = cos_pitch_angle(:u1, :B))
        @transform!(:w0 = clamp.(:w0, -1, 1), :w1 = clamp.(:w1, -1, 1))
        @transform!(:α0 = acosd.(:w0), :α1 = acosd.(:w1))
        @transform!(:Δα = :α1 .- :α0, :Δw = :w1 .- :w0)
        select!(Not(:wϕ0))
    end
end

function get_result(r::DataFrameRow, params=[:θ, :β, :v, :tmax, :B, :alg, :sign])
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