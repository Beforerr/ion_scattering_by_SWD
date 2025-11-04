using StaticArrays
using DrWatson
using DataFrames, DataFramesMeta
using CurrentSheetTestParticle: field_lines_asym_distance
# %% Solution processing

function extract_info(sol, B)
    u0 = SVector(sol.prob.u0...)
    u1 = SVector(sol.u[end]...)
    t1 = sol.t[end]
    dR_perp_asym = field_lines_asym_distance(sol)
    return (; u0, u1, t1, dR_perp_asym) |> pairs |> Dict
end

function process_sols(prob, sol, wϕs)
    results = extract_info.(sol.u, prob.B) |> DataFrame
    results.wϕ0 = wϕs
    results.B .= prob.B
    results.v .= prob.v
    results.tmax .= prob.tspan[2]
    return process_result!(results)
end

function process_result!(df::AbstractDataFrame)

    "wϕ0" in names(df) && @rtransform!(df, :μ0 = :wϕ0[1], :ϕ0 = :wϕ0[2])
    "dR_perp_asym" in names(df) && "v" in names(df) && @transform!(
        df,
        :dR_perp_asym_norm = :dR_perp_asym ./ :v
    )

    return @chain df begin
        @rtransform!(:μ1 = cos_pitch_angle(:u1, :B))
        @transform!(:μ0 = clamp.(:μ0, -1, 1), :μ1 = clamp.(:μ1, -1, 1))
        @transform!(:α0 = acosd.(:μ0), :α1 = acosd.(:μ1))
        @transform!(:s2α0 = sind.(:α0) .^ 2, :s2α1 = sind.(:α1) .^ 2)
        @transform!(:Δα = :α1 .- :α0, :Δμ = :μ1 .- :μ0, :Δs2α = :s2α1 - :s2α0)
        @transform!(:leave = :t1 .!= :tmax)
        select!(Not(:wϕ0))
    end
end
