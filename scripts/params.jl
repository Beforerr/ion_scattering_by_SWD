function example_params(; vs=[1, 8], θs=[60, 85], βs=[50, 75])
    allparams = Dict(
        :θ => θs,
        :β => βs,
        :v => vs,
        :tspan => (0, 1024),
        :init_kwargs => (; Nw=180, Nϕ=120),
    )
    return dict_list(allparams)
end

function test_params()
    θs = 45:20:85 |> collect
    ws = 60:60:180 |> collect
    βs = ws ./ 2
    vs = 8.0 .^ (0:3)
    diffeq = (; abstol=1e-5, reltol=1e-5, maxiters=1e6)

    allparams = Dict(
        :θ => θs,
        :β => βs,
        :v => vs,
        :init_kwargs => (; Nw=90, Nϕ=120),
        :tspan => (0, 1024),
        :diffeq => diffeq
    )

    return dict_list(allparams)
end

function test_sign_params()
    allparams = Dict(
        :θ => 65,
        :β => 60,
        :sign => [-1, 1],
        :v => 8.0,
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
        :diffeq => diffeq
    )

    return dict_list(allparams)
end

"""
Parameters for the simulation with tangential discontinuity

# Notes
This function is not updated to the latest version of the code.
"""
function test_params_TD()
    vs = 2.0 .^ (1:6)

    allparams = Dict(
        :Bfn => [TD_B_field, RD_B_field],
        :θ => 85,
        :β => 47.5,
        :v => vs,
        :init_kwargs => (; Nw=90, Nϕ=120),
    )

    return dict_list(allparams)
end