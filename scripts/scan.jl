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
    else
        return Vern9()
    end
end

function makesim(d::Dict; save_everystep = false, kwargs...)
    @unpack θ, β, v, alg_sym, init_kwargs, sign, tspan = d
    B = RD_B_field(; θ, β, sign)
    wϕs = w_ϕ_pairs(; init_kwargs...)
    filter_wϕs!(wϕs, θ)
    u0s = init_state(B, v, wϕs)

    isoutofdomain = isoutofdomain_params(v)
    alg = sym2ins(alg_sym)
    sol = solve_params(B, u0s; alg, tspan, save_everystep, isoutofdomain)
    results = extract_info.(sol.u) |> DataFrame
    results.wϕ0 = wϕs
    return @dict results α β sign v alg_sym
end

function main()
    θs = deg2rad.(5 : 10 : 175) # from 5° to 175° in 10° steps
    ws = 25:10:175 # PDF is reliable at β > 15°, corresponding to rotation angle $w$ from 30° to 180° 
    βs = deg2rad.(ws ./ 2)
    vs = 2. .^ (-2:8)

    allparams = Dict(
        :θ => θs,
        :β => βs,
        :sign => [-1, 1],
        :v => vs,
        :alg_sym => [:AutoVern9],
        :init_kwargs => (;Nw=90, Nϕ=120),
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

function test_save_efficiency(;test_dict = first(dicts))
    path = datadir("simulations")
    @info "Running simulation" test_dict
    result, file = produce_or_load(makesim, test_dict, path)
    # Check the filesize
    filesize = stat(file).size
    @info "Filesize" filesize
    return result, file
end
