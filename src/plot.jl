using AlgebraOfGraphics
using CairoMakie
using DataFrames, DataFramesMeta
using Match

begin
    set_aog_theme!()
    theme = (; colormap=(:batlow))
    update_theme!(; theme...)
end

vals(df, s) = unique(df[!, s]) |> sort

begin
    w0 = :w0 => "cos(α₀)"
    w1 = :w1 => "cos(α₁)"
    Δw = :Δw => "Δ cos α"
    α0 = :α0 => "α₀"
    α1 = :α1 => "α₁"
    Δα = :Δα
    ϕ0 = :ϕ0 => "ϕ₀"
    xyα = (α0, α1)
    xyw = (w0, w1)
    Δt = :t1 => "Δt"
    leave = :leave => renamer([true => "Leave", false => "Trapped"])

    ergs_approx = ["~10 eV", "~100 eV", "~5 keV", "~100 keV", "~1 MeV"]
end

begin
    rename_sym(s, v) = "$(s) = $(round(Int, v))"
    rename_sym(s) = v -> rename_sym(s, v)
    rename_sym_deg(s) = v -> (rename_sym(s, v) * "°")
    get_map(s; renamer=rename_sym) = s => renamer(s)
    sign_map = :sign => renamer([1 => "Left-handed", -1 => "Right-handed"])
    θ_map = get_map(:θ; renamer=rename_sym_deg)
    β_map = get_map(:β; renamer=rename_sym_deg)
    v_map = get_map(:v)
end


begin
    colorscale = log10
    # density_layer = AlgebraOfGraphics.density(npoints=32)
    density_layer() = histogram(; bins=64, normalization=:pdf) * visual(; colorscale)
    colorrange(; scale=colorscale) = scale == log10 ? (1e-2, 1e1) : (0, 5)
    tm_scale() = scales(Color=(;
        colorrange=colorrange()
    ))
end

w_axis = (; limits=((-1, 1), (-1, 1)), yreversed=true, yticks=[-1, 0, 1], xticks=[-1, 0, 1])
α_axis = (; limits=((0, 180), (0, 180)))

function pa_pair_plot(layer; axis=(;), kwargs...)
    v = visual(alpha=0.1, markersize=3; legend=(; markersize=12))
    l = layer * mapping(w0, w1) * v
    draw(l; axis, kwargs...)
end

pa_pair_plot(df::AbstractDataFrame) = pa_pair_plot(data(df))

# Function to plot one figure per column
function pa_layer(s)
    layout = @match s begin
        :v => mapping(col=θ_map, row=β_map)
        :β => mapping(col=θ_map, row=v_map)
    end
    return mapping(xyw...) * density_layer() * layout
end

function pa_pair_hist!(df::AbstractDataFrame, layout, s=:v; layer=pa_layer(s), scales=tm_scale(), axis=w_axis, kwargs...)
    vs = vals(df, s)  # Get unique velocity values
    fgs = [layout[1, i] for i in 1:length(vs)]
    return map(zip(fgs, vs)) do (fg, v)
        df_s = @subset(df, $s .== v)
        plt = data(df_s) * layer
        grids = draw!(fg, plt, scales; axis, kwargs...)
        r = s == :v ? rename_sym : rename_sym_deg
        Label(fg[0, :], r(s)(v), tellwidth=false)
        grids
    end
end

function pa_pair_hist(df; figure=(; size=(1200, 400)), scale=colorscale, kw...)
    fig = Figure(; figure...)
    grids = pa_pair_hist!(df, fig[1, 1]; kw...)
    colorbar!(fig[1:end, end+1], grids[1]; scale)
    fig, grids
end

function sign_label(fig; i1=0, i2=2)
    Label(fig[i1, :], "Left-hand rotation", font=:bold, tellwidth=false)
    Label(fig[i2, :], "Right-hand rotation", font=:bold, tellwidth=false)
end

function pa_pair_hist(ldf, rdf; figure=(; size=(1200, 800)), scale=colorscale, kw...)
    fig = Figure(; figure...)
    grids = pa_pair_hist!(ldf, fig[1, 1]; kw...)
    grids = pa_pair_hist!(rdf, fig[3, 1]; kw...)
    sign_label(fig)
    colorbar!(fig[1:end, end+1], grids[1]; scale)
    fig
end

# %%
"""
Plot the distribution of pitch-angle cosine variation
"""
function pa_diff_plot!(fig, layer; kwargs...)
    v = AlgebraOfGraphics.density() * visual(Lines)
    plt = layer * mapping(Δα) * v
    axis = (; xlabel="Δα", ylabel="f (Δα)", limits=((-2, 1), (0, 4)))
    draw!(fig, plt; axis, kwargs...)
end

function pa_diff_plot(layer; kwargs...)
    v = AlgebraOfGraphics.density() * visual(Lines)
    plt = layer * mapping(Δα) * v

    axis = (; xlabel="Δα", ylabel="f (Δα)", limits=((-2, 1), (0, 4)))
    draw(plt; axis, kwargs...)
end

pa_diff_plot!(fig, df::AbstractDataFrame; kwargs...) = pa_diff_plot!(fig, data(df) * mapping(col=θ_map, row=β_map, color=v_map); kwargs...)
pa_diff_plot(df::AbstractDataFrame; kwargs...) = pa_diff_plot(data(df) * mapping(col=θ_map, row=β_map, color=v_map); kwargs...)


# Create Two-dimensional maps of the final value w1 are plotted as functions of the initial pitch-angle cosine w0 and gyrophase φ0
function w1_map_plot(l; color=w1, scale=scales(;), kwargs...)
    plt = l * mapping(w0, ϕ0, color=color)
    draw(plt, scale; kwargs...)
end

#%%

function get_ax(layout, idxs)
    is3D = length(idxs) >= 3
    pos = layout[1, 1]
    ax = is3D ? Axis3(pos) : Axis(pos)
    is3D && (ax.protrusions = 50) # removes overlap of labels
    return ax
end

function create_slider(parameters, layout)
    tuples_for_slidergrid = map(collect(parameters)) do (label, range)
        label = string(label)
        (; label=label, range=range, startvalue=first(range))
    end
    sg = SliderGrid(layout[1, 1], tuples_for_slidergrid...)
    return sg.sliders
end

function plot_params(parameter_sliders, idxs)
    fig = Figure()
    layout = fig[1, 1] = GridLayout()
    paramlayout = fig[2, :] = GridLayout(tellheight=true, tellwidth=false)

    sliders = create_slider(parameter_sliders, paramlayout)
    sliders_values = map(s -> s.value, sliders)

    ax = get_ax(layout, idxs)

    lift(sliders_values...) do vs...
        empty!(ax)
        sols = solve_params(α, β, vs...)
        for sol in sols.u
            lines!(ax, sol[idxs, :])
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
        plot!(sol_field, idxs=(1, 2, 3), label="Magnetic Field Line", color=:gray, alpha=0.5)
    end

    return f
end

E0(x) = SVector(0.0, 0.0, 0.0)

function get_gc_func(B)
    param = prepare(E0, B, species=User)
    get_gc(param)
end

function plot_gc!(sol, B; color=Makie.wong_colors()[2])
    gc = get_gc_func(B)
    gc_plot(x, y, z, vx, vy, vz) = (gc([x, y, z, vx, vy, vz])...,)
    lines!(sol, idxs=(gc_plot, 1, 2, 3, 4, 5, 6); color)
end

function plot_gc_field_lines!(sol, B; idxs=(1, 2, 3), kwargs...)
    gc = get_gc_func(B)
    gc0 = gc(sol[1])
    gcf = gc(sol[end])

    isoutofdomain = (u, p, t) -> abs(u[3]) > maximum(abs.(sol[3, :]))

    tmax = sol.t[end]
    fl0_sol = CurrentSheetTestParticle.solve_fl(gc0, B; tspan=(0.0, tmax), isoutofdomain, kwargs...)
    flf_sol = CurrentSheetTestParticle.solve_fl(gcf, B; tspan=(0.0, -tmax), isoutofdomain, kwargs...)
    lines!(fl0_sol; idxs, color=Makie.wong_colors()[3])
    lines!(flf_sol; idxs, color=Makie.wong_colors()[4])
end