using Revise
using CurrentSheetTestParticle
using DrWatson
using DataFrames
using Logging, LoggingExtras
using ProgressMeter
using StaticArrays
include("./params.jl")
include("../src/utils.jl")

filtered_logger = ActiveFilteredLogger(global_logger()) do args
    # filter "you passed a key as a symbol instead of a string"
    re_JLD2 = r"you passed a key"
    re_SCIML = r"dt was forced below floating point"
    return !occursin(re_JLD2, args.message) && !occursin(re_SCIML, args.message)
end

error_only_logger = MinLevelLogger(current_logger(), Logging.Error);

function makesim(d; save_everystep=false, kwargs...)
    @unpack θ, β, v, sign, alg, init_kwargs, diffeq, tspan = ProblemParams(; d...)
    B = RD_B_field(; θ, β, sign)
    u0s, wϕs = init_states_pm(B, v; init_kwargs...)

    isoutofdomain = isoutofdomain_params(v)
    sol = solve_params(B, u0s; alg, tspan, diffeq, save_everystep, isoutofdomain)
    results = extract_info.(sol.u) |> DataFrame
    results.wϕ0 = wϕs
    return merge(d, @dict results)
end

"""
Parameters with high resolution
"""
function scan_params_hi()
    vs = 2.0 .^ (1:6)

    allparams = Dict(
        :θ => 85,
        :β => 47.5,
        :v => vs,
        :tspan => (0, 1024),
        :init_kwargs => (; Nw=180, Nϕ=120),
    )
    return dict_list(allparams)
end

function scan_params()
    θs = 5:10:85 |> collect # from 5° to 85 in 10° steps
    ws = 25:10:175 |> collect # PDF is reliable at β > 15°, corresponding to rotation angle $w$ from 30° to 180° 
    βs = ws ./ 2
    vs = 2.0 .^ (-2:8)

    allparams = Dict(
        :θ => θs,
        :β => βs,
        :v => vs,
        :init_kwargs => (; Nw=90, Nϕ=120),
        :tspan => (0, 1024),
    )
    return dict_list(allparams)
end

function main(; params=scan_params(), dir="simulations", logger=error_only_logger)
    path = datadir(dir)
    with_logger(logger) do
        @showprogress map(params) do d
            produce_or_load(makesim, d, path; loadfile=false)
        end
    end
end

test(; params=test_params(), dir="test", logger=filtered_logger) = main(; params, dir, logger)
test_alg(; params=test_params_alg(), dir="test_alg", logger=filtered_logger) = main(; params, dir, logger)

(@main)(ARGS) = main()

function test_save_efficiency(; test_dict=first(dicts))
    path = datadir("simulations")
    @info "Running simulation" test_dict
    result, file = produce_or_load(makesim, test_dict, path)
    # Check the filesize
    filesize = stat(file).size
    @info "Filesize" filesize
    return result, file
end
