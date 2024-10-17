using Revise
using CurrentSheetTestParticle
using DrWatson
using DataFrames
using Logging, LoggingExtras
using ProgressMeter
using StaticArrays
using OrdinaryDiffEq
include("./params.jl")

filtered_logger = ActiveFilteredLogger(global_logger()) do args
    # filter "you passed a key as a symbol instead of a string"
    re_JLD2 = r"you passed a key"
    re_SCIML = r"dt was forced below floating point"
    return !occursin(re_JLD2, args.message) && !occursin(re_SCIML, args.message)
end

error_only_logger = MinLevelLogger(current_logger(), Logging.Error);

function extract_info(sol)
    u0 = SVector(sol.prob.u0...)
    u1 = SVector(sol.u[end]...)
    t1 = sol.t[end]
    return @dict u0 u1 t1
end

function sym2ins(alg)
    if alg == :AutoVern9
        return AutoVern9(Rodas4P())
    elseif alg == :Boris
        @warn "Boris algorithm is not implemented yet"
    else
        return Vern9()
    end
end

function makesim(d::Dict; save_everystep=false, kwargs...)
    @unpack θ, β, v, sign, alg_sym, init_kwargs, diffeq, tspan = d
    θ, β = deg2rad.([θ, β])
    B = RD_B_field(; θ, β, sign)
    u0s, wϕs = init_states_pm(B, v; init_kwargs...)

    isoutofdomain = isoutofdomain_params(v)
    alg = sym2ins(alg_sym)
    sol = solve_params(B, u0s; alg, tspan, diffeq, save_everystep, isoutofdomain)
    results = extract_info.(sol.u) |> DataFrame
    results.wϕ0 = wϕs
    return merge(d, @dict results)
end

function scan_params()
    θs = 5:10:85 # from 5° to 85 in 10° steps
    ws = 25:10:175 # PDF is reliable at β > 15°, corresponding to rotation angle $w$ from 30° to 180° 
    βs = ws ./ 2
    vs = 2.0 .^ (-2:8)
    diffeq = CurrentSheetTestParticle.DEFAULT_DIFFEQ_KWARGS

    allparams = Dict(
        :θ => θs,
        :β => βs,
        :sign => [1],
        :v => vs,
        :alg_sym => [:AutoVern9],
        :init_kwargs => (; Nw=90, Nϕ=120),
        :tspan => (0, 1024),
        :diffeq => diffeq
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
