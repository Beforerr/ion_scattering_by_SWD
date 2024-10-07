using CairoMakie

function plot_trajectory(sol, sol_field=missing, sol_gc=missing)
    # Plot the actual trajectory
    f = Figure()

    ax = Axis3(
        f[1, 1],
        title="Test Particle Trajectory",
        xlabel="x",
        ylabel="y",
        zlabel="z",
        aspect=:data,
    )

    plot!(sol, idxs=(1, 2, 3))

    # Plot the guiding center trajectory if computed
    if sol_gc !== missing
        plot!(sol_gc, idxs=(1, 2, 3), color=:red, label="Guiding Center")
    end

    # Plot Magnetic Field Line
    if sol_field !== missing
        plot!(sol_field, idxs=(1, 2, 3), label="Magnetic Field Line", color = :gray, alpha = 0.5)
    end

    return f
end