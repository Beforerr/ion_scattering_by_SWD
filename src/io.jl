using DataFramesMeta
using CurrentSheetTestParticle

#%%
subset_leave(df) = @subset(df, :t1 .!= :tmax)
subset_trap(df) = @subset(df, :t1 .== :tmax)

function get_result_dfs(; dir="simulations")
    @chain begin
        collect_results(datadir(dir))
        @rtransform(
            :B = RD_B_field(; θ=:θ, β=:β, sign=:sign)
        )
        @transform!(
            :θ = round.(rad2deg.(:θ)), 
            :β = round.(rad2deg.(:β)),
            :tmax = last.(:tspan)
        )
        select!(Not(:tspan))
    end 
end

function get_result(; dir="simulations")
    dfs = get_result_dfs(; dir)
    results = map(get_result, eachrow(dfs))
    vcat(results...)
end

function get_result(r::DataFrameRow)
    params = [:θ, :β, :sign, :v, :alg_sym, :tmax, :B]
    @chain df begin
        insertcols!(
            Pair.(params, Array(r[params]))...; makeunique=true
        )
        @rtransform!(:w0 = :wϕ0[1], :ϕ0 = :wϕ0[2], :w1 = cos_pitch_angle(:u1, :B))
        @transform!(:α0 = acosd.(:w0), :α1 = acosd.(:w1))
        @transform!(:Δα = :α1 .- :α0)
        select!(Not(:wϕ0))
    end
end

function split_results(df)
    ldf = @subset(df, :sign .== 1)
    rdf = @subset(df, :sign .== -1)
    return ldf, rdf
end