using DrWatson
using AlgebraOfGraphics
using DataFrames
using DataFramesMeta
using CurrentSheetTestParticle
using CairoMakie

set_aog_theme!()
theme = (;colormap = Reverse(:viridis))
update_theme!(; theme...)

params = [:α, :β, :v, :sign]

function get_results()
    dfs = collect_results(datadir("simulations"))

    for r in eachrow(dfs)
        α, β, sign = r[:α], r[:β], r[:sign]
        insertcols!(r[:results], Pair.(params, Array(r[params]))...; makeunique=true)
        insertcols!(r[:results], :B => RD_B_field(; α, β, sign))
    end

    @chain begin
        vcat(dfs[!, :results]...)
        @rename!(:θ = :α)
        @rtransform!(:w0 = :wϕ0[1], :ϕ0 = :wϕ0[2], :w1 = cos_pitch_angle(:u1, :B))
        @transform!(:Δα = :w1 .- :w0,)
        transform!([:θ, :β] .=> ByRow(rad2deg), renamecols=false)
        select!(Not(:wϕ0))
    end
end

function split_results(df)
    ldf = @subset(df, :sign .== 1)
    rdf = @subset(df, :sign .== -1)
    return ldf, rdf
end

subset_v(df) = @subset(df, :v .<= 128)
subset_v(v::Number) = df -> @subset(df, :v .<= v)
subset_β(itr) = df -> @rsubset(df, :β ∈ itr)

#%%
const DEFAULT_TMAX = CurrentSheetTestParticle.DEFAULT_TSPAN[2]
subset_leave(df; tmax=DEFAULT_TMAX) = @subset(df, :t1 .!= tmax)

vals(df, s) = unique(df[!, s]) |> sort
rename_sym(s, v) = "$(s) = $(round(Int, v))"
rename_sym(s) = v -> rename_sym(s, v)
rename_sym_deg(s) = v -> (rename_sym(s, v) * "°")
get_map(s; renamer = rename_sym) = s => renamer(s)

θ_map = get_map(:θ; renamer = rename_sym_deg)
β_map = get_map(:β; renamer = rename_sym_deg)
v_map = get_map(:v)
sign_map = :sign => renamer([1 => "Left-handed", -1 => "Right-handed"])
Δα = :Δα
w0 = :w0 => "α₀"
w1 = :w1 => "α₁"
ϕ0 = :ϕ0 => "ϕ₀"

#%%
function pa_pair_plot(layer; axis=(;), kwargs...)
    v = visual(alpha=0.1, markersize=3; legend=(; markersize=12))
    l = layer * mapping(w0, w1) * v
    draw(l; axis, kwargs...)
end

w_axis = (; limits=((0, 1), (-1, 1)))

pa_pair_plot(df::AbstractDataFrame) = pa_pair_plot(data(df))

# Function to plot one figure per column
function pa_pair_plot_by_col!(df::AbstractDataFrame, s, axs, layer; scale=scales(), axis=w_axis, kwargs...)
    vs = (sort ∘ unique)(df[!, s])  # Get unique velocity values
    return map(zip(axs, vs)) do (fg, v)
        df_s = @subset(df, $s .== v)
        plt = data(df_s) * mapping(w0, w1) * layer
        grids = draw!(fg, plt, scale; axis, kwargs...)
        r =  s == :v ? rename_sym : rename_sym_deg
        Label(fg[0, :], r(s)(v))
        grids
    end
end

pa_pair_plot_by_v!(df, axs) = pa_pair_plot_by_col!(df, :v, axs, mapping(col=θ_map, row=β_map) * density_layer; scale)
pa_pair_plot_by_β!(df, axs) = pa_pair_plot_by_col!(df, :β, axs, mapping(col=θ_map, row=v_map) * density_layer; scale)

colorbar!

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

function sign_label(fig; i1=0, i2=2)
    Label(fig[i1, :], "Left-hand rotation", font = :bold, tellwidth=false)
    Label(fig[i2, :], "Right-hand rotation", font = :bold, tellwidth=false)
end

begin
    # density_layer = AlgebraOfGraphics.density(npoints=32)
    density_layer = histogram(; bins=64, normalization=:pdf)
    scale = scales(Color=(;
        nan_color=:transparent,
        lowclip=:transparent,
        colorrange=(0.001, 10)
    ))
end


# Define the colors for the colormap
using ColorSchemes

cmap0 = :BrBg
cs = colorschemes[cmap0]
Δw_cmap = cgrad([cs[1], get(cs, 0.5), get(cs, 0.75)], [0.0, 0.67, 1.0])
Δw_scale = scales(Color=(; colormap=Δw_cmap))

# Create Two-dimensional maps of the final value w1 are plotted as functions of the initial pitch-angle cosine w0 and gyrophase φ0, for α = 45◦ and β = 15◦ (upper panel), β = 45◦ (middle panel), β = 90◦ (lower panel).
function w1_map_plot(l; color = w1, scale= scales(;), kwargs...)
    plt = l * mapping(w0, ϕ0, color=color)
    draw(plt, scale; kwargs...)
end