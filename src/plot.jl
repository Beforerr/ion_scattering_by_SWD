function get_ax(layout, idxs)
    is3D = length(idxs) == 3
    pos = layout[1,1]
    ax = is3D ?  Axis3(pos) : Axis(pos)
    is3D && (ax.protrusions = 50) # removes overlap of labels
    return ax
end

function create_slider(parameters, layout)
    tuples_for_slidergrid = map(collect(parameters)) do (label, range)
        label = string(label)
        (; label=label, range = range, startvalue = first(range))
    end
    sg = SliderGrid(layout[1,1], tuples_for_slidergrid...)
    return sg.sliders
end

function plot_params(parameter_sliders, idxs)
    fig = Figure()
    layout = fig[1,1] = GridLayout()
    paramlayout = fig[2, :] = GridLayout(tellheight = true, tellwidth = false)

    sliders = create_slider(parameter_sliders, paramlayout)
    sliders_values = map(s -> s.value, sliders)

    ax = get_ax(layout, idxs)
    
    lift(sliders_values...) do vs...
        empty!(ax)
        sols = solve_params(α, β, vs...)
        for sol in sols.u
            lines!(ax,sol[idxs, :])
        end
    end
    fig
end

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