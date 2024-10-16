using Revise
using CurrentSheetTestParticle
using DrWatson
using DataFrames
using Logging, LoggingExtras
using ProgressMeter
using StaticArrays
using DifferentialEquations

logger = ActiveFilteredLogger(global_logger()) do args
    # filter "you passed a key as a symbol instead of a string"
    re_JLD2 = r"you passed a key"
    re_SCIML = r"dt was forced below floating point"
    return !occursin(re_JLD2, args.message) && !occursin(re_SCIML, args.message)
end

svec(x) = SVector(x...)

function extract_info(sol)
    u0 = sol.prob.u0 |> svec
    t = sol.t |> svec
    u = map(svec, sol.u)

    u1 = u[end]
    t1 = t[end]
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
    @unpack θ, β, v, sign, alg_sym, init_kwargs, tspan = d
    B = RD_B_field(; θ, β, sign)
    wϕs = w_ϕ_pairs(; init_kwargs...)
    filter_wϕs!(wϕs, θ)
    u0s = init_state(B, v, wϕs)

    isoutofdomain = isoutofdomain_params(v)
    alg = sym2ins(alg_sym)
    sol = solve_params(B, u0s; alg, tspan, save_everystep, isoutofdomain)
    results = extract_info.(sol.u) |> DataFrame
    results.wϕ0 = wϕs
    return merge(d, @dict results)
end

function test_params()
    θs = deg2rad.(50:20:130) # from 50° to 130° in 20° steps
    ws = 55:30:175
    βs = deg2rad.(ws ./ 2)
    vs = 4.0 .^ (0:4)

    allparams = Dict(
        :θ => θs,
        :β => βs,
        :sign => [-1, 1],
        :v => vs,
        :alg_sym => [:AutoVern9, :Boris],
        :init_kwargs => (; Nw=90, Nϕ=120),
        :tspan => (0, 256),
    )

    return dict_list(allparams)
end

function test()
    dicts = test_params()

    @showprogress map(dicts) do d
        path = datadir("test")
        with_logger(logger) do
            produce_or_load(makesim, d, path; loadfile=false)
        end
    end
end

function main()
    θs = deg2rad.(5:10:175) # from 5° to 175° in 10° steps
    ws = 25:10:175 # PDF is reliable at β > 15°, corresponding to rotation angle $w$ from 30° to 180° 
    βs = deg2rad.(ws ./ 2)
    vs = 2.0 .^ (-2:8)

    allparams = Dict(
        :θ => θs,
        :β => βs,
        :sign => [-1, 1],
        :v => vs,
        :alg_sym => [:AutoVern9],
        :init_kwargs => (; Nw=90, Nϕ=120),
        :tspan => (0, 512),
    )

    dicts = dict_list(allparams)

    @showprogress map(dicts) do d
        path = datadir("simulations")
        with_logger(logger) do
            produce_or_load(makesim, d, path; loadfile=false)
        end
    end
end

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
