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