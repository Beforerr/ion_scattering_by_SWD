function test_params()
    θs = 45:20:85
    ws = 60:60:180
    βs = ws ./ 2
    vs = 8.0 .^ (0:2)
    diffeq = (; abstol=1e-5, reltol=1e-5, maxiters=1e6)

    allparams = Dict(
        :θ => θs,
        :β => βs,
        :sign => [-1, 1],
        :v => vs,
        :alg_sym => [:AutoVern9],
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
    θs = 45:20:85
    ws = 60:60:180
    βs = ws ./ 2
    vs = 8.0 .^ (0:2)
    diffeq = (; abstol=1e-5, reltol=1e-5, maxiters=1e5)

    allparams = Dict(
        :θ => θs,
        :β => βs,
        :sign => [-1, 1],
        :v => vs,
        :alg_sym => [:ImplicitMidpoint, :AutoVern9],
        :init_kwargs => (; Nw=90, Nϕ=120),
        :tspan => (0, 512),
        :diffeq => diffeq
    )

    return dict_list(allparams)
end