function test_params()
    θs = 45:20:85 |> collect
    ws = 60:60:180 |> collect
    βs = ws ./ 2
    vs = 8.0 .^ (0:2)
    diffeq = (; abstol=1e-5, reltol=1e-5, maxiters=1e6)

    allparams = Dict(
        :θ => θs,
        :β => βs,
        :sign => [-1, 1],
        :v => vs,
        :init_kwargs => (; Nw=90, Nϕ=120),
        :tspan => (0, 1024),
        :diffeq => diffeq
    )

    return dict_list(allparams)
end

"""
Test parameters for the simulation with different algorithms
"""
function test_params_alg()
    θs = 45:20:85 |> collect
    ws = 60:60:180 |> collect
    βs = ws ./ 2
    vs = 8.0 .^ (0:2)
    diffeq = (; abstol=1e-5, reltol=1e-5, maxiters=1e5)

    allparams = Dict(
        :θ => θs,
        :β => βs,
        :v => vs,
        :alg => [:ImplicitMidpoint, :AutoVern9],
        :init_kwargs => (; Nw=90, Nϕ=120),
        :tspan => (0, 512),
        :diffeq => diffeq
    )

    return dict_list(allparams)
end