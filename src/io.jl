using DataFramesMeta

#%%
subset_leave(df) = @subset(df, :t1 .!= :tmax)
subset_trap(df) = @subset(df, :t1 .== :tmax)

function get_results(; dir="simulations")
    dfs = collect_results(datadir(dir))

    for r in eachrow(dfs)
        get_result(r)
    end

    @chain begin
        vcat(dfs[!, :results]...)
        @transform!(:α0 = acosd.(:w0), :α1 = acosd.(:w1))
        @transform!(:Δα = :α1 .- :α0)
    end
end

function get_result(r)
    @unpack θ, β, sign, v, alg_sym, tspan = r
    params = [:θ, :β, :sign, :v, :alg_sym]
    df = copy(r[:results])
    @chain df begin
        insertcols!(
            :B => RD_B_field(; θ, β, sign),
            :tmax => tspan[2],
            Pair.(params, Array(r[params]))...; makeunique=true
        )
        @rtransform!(:w0 = :wϕ0[1], :ϕ0 = :wϕ0[2], :w1 = cos_pitch_angle(:u1, :B))
        transform!([:θ, :β] .=> ByRow(rad2deg), renamecols=false)
        select!(Not(:wϕ0))
    end
end

function split_results(df)
    ldf = @subset(df, :sign .== 1)
    rdf = @subset(df, :sign .== -1)
    return ldf, rdf
end