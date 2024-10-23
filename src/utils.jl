using StaticArrays
using DrWatson
using DataFrames, DataFramesMeta

# %% Solution processing

function extract_info(sol)
    u0 = SVector(sol.prob.u0...)
    u1 = SVector(sol.u[end]...)
    t1 = sol.t[end]
    return @dict u0 u1 t1
end

function process_sols(sol, B, wϕs)
    results = extract_info.(sol.u) |> DataFrame
    results.wϕ0 = wϕs
    results.B .= B
    process_result!(results)
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