using Revise
includet("../src/CurrentSheetTestParticle.jl")
using DrWatson
using LinearAlgebra
using OrdinaryDiffEq
using DataFrames
using Logging
using ProgressMeter
using StaticArrays

svec(x) = SVector(x...)

function extract_info(sol)
    u0 = sol.prob.u0 |> svec
    t = sol.t |> svec
    u = map(svec, sol.u)

    u1 = u[end]
    t1 = t[end]
    return @dict u0 u1 t1
end

function makesim(d::Dict; save_everystep = false, kwargs...)
    @unpack α, β, v, alg, init_kwargs, sign = d
    alg = eval(alg)()
    B = r -> RD_B_field(r, α, β; sign)
    wϕs = w_ϕ_pairs(; init_kwargs...)
    u0s = init_state(B, v, wϕs)

    sol = solve_params(B, v, u0s; alg, save_everystep)
    results = extract_info.(sol.u) |> DataFrame
    results.wϕ0 = wϕs
    return @dict results α β sign v alg
end

Logging.disable_logging(Logging.Warn)

function main()
    dα = π / 16
    αs = collect(π / 4 : dα : π / 2 - dα)
    βs = collect(π / 12 : π / 12 : π / 2)
    vs = [1, 8, 64, 128, 1024, 8192]
    allparams = Dict(
        :α => αs,
        :β => βs,
        :sign => [-1, 1],
        :v => vs,
        :alg => [:Vern9],
        :init_kwargs => (;Nw=64, Nϕ=128),
        :tspan => (0, 256),
    )
       
    dicts = dict_list(allparams)

    @showprogress map(dicts) do d
        path = "./data/simulations"
        @info "Running simulation" d
        produce_or_load(makesim, d, path; loadfile=false)
    end
end

(@main)(ARGS) = main()

function test_save_efficiency(;test_dict = first(dicts))
    path = "./data/simulations"
    @info "Running simulation" test_dict
    result, file = produce_or_load(makesim, test_dict, path)
    # Check the filesize
    filesize = stat(file).size
    @info "Filesize" filesize
    return result, file
end
