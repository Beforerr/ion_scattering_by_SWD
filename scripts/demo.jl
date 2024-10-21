using Revise
using DrWatson
using CurrentSheetTestParticle
using Beforerr
include("../src/plot.jl")
include("../src/utils.jl")

# ---------------------------
# Minimal working example
# ---------------------------
d = ProblemParams(
    θ=45,
    β=90,
    v=1,
    init_kwargs = (; w=[-0.8,0.8])
)

sols, (wϕs, B) = solve_params(d);

Bx(z) = B([0, 0, z])[1]
By(z) = B([0, 0, z])[2]
Bz(z) = B([0, 0, z])[3]

using GLMakie
GLMakie.activate!()

θ = d.θ
ku(x; θ=θ) = cosd(θ) * x
ku(t, x) = (t, ku(x))
ku(x, y, z) = (ku(x), ku(y), z)
ku_zx(z, x) = (z, ku(x))

Eₖ(u) = 1/2 * sum(u[4:6].^2)

function plot_sols(sols, idxs)
    fig = Figure()
    layout = fig[1, 1]
    ax = get_ax(layout, idxs)
    for sol in sols.u
        plot!(ax, sol, idxs=idxs)
    end
    fig
end

# Plot trajectories of three particles projected onto the z-x plane.
function plot_detail(sols; idxs = (3, 2))
    fig = Figure(; size=(800, 400))
    ax1 = Axis(fig[1, 1]; xlabel="z", ylabel="B", limits = (nothing, (-1, 1)))
    ax2 = Axis(fig[2, 1]; xlabel="z", ylabel="y")
    ax3 = Axis(fig[3, 1]; xlabel="z", ylabel="x")
    # ax3 = Axis(fig[3, 1]; xlabel="x", ylabel="y")
    ax4 = Axis(fig[1:end, 2]; xlabel = "t", ylabel = "Cos(α)")
    # ax5 = Axis(fig[2, 2]; limits = (nothing, (0, 0.5)))

    zmin, zmax = minimum(sols[1][3, :]), maximum(sols[1][3, :])
    z = range(zmin, zmax, length=100)
    lines!(ax1, z, Bx, label="Bx")
    lines!(ax1, z, By, label="By")
    lines!(ax1, z, Bz, label="Bz")
    axislegend(ax1, orientation = :horizontal, position = :rb)

    for sol in sols.u
        # plot!(ax2, sol, idxs=(3, 2)) # this would smooth the trajectory
        lines!(ax2, sol[3, :], sol[2, :])
        lines!(ax3, sol[3, :], sol[1, :])
        # plot!(ax3, sol, idxs=(1,2))
        scatterlines!(ax4, sol.t, pa.(sol.u, B), alpha=0.5)
        # scatterlines!(ax5, sol.t, Eₖ.(sol.u), alpha=0.5)
    end
    fig
end

idxs = (ku, 1, 2, 3)
plot_sols(sols[16:18], idxs)

plot_detail(sols[[12,13]])
easy_save("example_tp", plot_detail(sols[[12,13]]))

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

# Histogram
d = ProblemParams(
    θ=85,
    β=42.5,
    v=2,
    init_kwargs=(; Nw=90, Nϕ=120)
)

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