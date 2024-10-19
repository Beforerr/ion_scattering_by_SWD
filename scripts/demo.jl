using Revise
using DrWatson
using CurrentSheetTestParticle
include("../src/plot.jl")
include("../src/utils.jl")

# ---------------------------
# Minimal working example
# ---------------------------
d = ProblemParams(
    θ=85,
    β=47.5,
    v=2,
    init_kwargs=(; Nw=90, Nϕ=120)
)

sol, (wϕs, B) = solve_params(d);

using GLMakie
GLMakie.activate!()

const pa = cos_pitch_angle

ku(x) = cos(θ) * x
ku(t, x) = (t, ku(x))
ku(x, y, z) = (cos(θ) * x, cos(θ) * y, z)

function plot_sols(sols, idxs)
    fig = Figure()
    layout = fig[1, 1]
    ax = get_ax(layout, idxs)
    for sol in sols.u
        plot!(ax, sol, idxs=idxs)
    end
    fig
end
idxs = (ku, 1, 2, 3)
plot_sols(sols, idxs)

observables = [3, 6, Eₖ, pa]

let sols = sols.u[1:8:end]
    fig = Figure()
    ax = Axis(fig[1, 1])
    for sol in sols
        scatterlines!(sol.t, pa.(sol.u, B), alpha=0.5)
    end
    ax = Axis(fig[2, 1])
    for sol in sols
        scatterlines!(sol.t, sol[3, :], alpha=0.5)
    end
    fig
end

let
    results = extract_info.(sol.u) |> DataFrame
    results.wϕ0 = wϕs
    results.B .= B
    process_result!(results)

    v = histogram(; bins=64)
    plt0 = data(results) * mapping(xyw...) * v
    fg = draw(plt0; axis=w_axis)
    Label(fg.figure[0, 1:end], "θ = $(d.θ)°, β = $(d.β)°")
    fg
end

# ---------------------------
α, β, ϕ = π / 4, π / 2, 0.125

vrange = logrange(64, 0.125, 8)
wrange = 0:0.1:1
v, w = 8, 0.5

parameter_sliders = OrderedDict(:v => vrange, :w => wrange)
idxs = [1, 2, 3]

sols = solve_params(α, β, v, w)
plot_params(parameter_sliders, idxs)

# ---------------------------
# Dynamical system
# ---------------------------
using DynamicalSystems

param = prepare(E_field, B_field; species=User)

u0s = init_state(α, β, v, w; Nϕ=64)

isoutofdomain = isoutofdomain_params(v)
diffeq = (; alg=alg, isoutofdomain)

total_time = 100
Y, t = trajectory(ds, total_time)


fig = Figure()
ax = Axis(fig[1, 1]; xlabel="time", ylabel="variable")
for var in columns(Y)
    lines!(ax, t, var)
end
fig

# compute the Lyapunov spectrum
steps = 10_000
lyapunovspectrum(ds, steps)
# As expected, there is at least one positive Lyapunov exponent, because the system is chaotic, and at least one zero Lyapunov exponent, because the system is continuous time.
z_max = 2 * init_z_pos(v) |> abs
timeseries_ylims = [(-z_max, z_max), (-v, v), missing, (-1, 1)]

lims = ((-z_max, z_max), (-v, v))
idxs = [3, 6]

fig, dsobs = interactive_trajectory_timeseries(
    ds, observables, u0s;
    lims, idxs,
    timeseries_ylims
)
fig


figure, oddata = interactive_orbitdiagram(ds, 3)