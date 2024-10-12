using DrWatson
using AlgebraOfGraphics
using CairoMakie
using DataFrames
using DataFramesMeta
using Revise
include("../src/CurrentSheetTestParticle.jl")
set_aog_theme!()
theme = (;colormap = :deep)
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
        @rtransform!(:w0 = :wϕ0[1], :ϕ0 = :wϕ0[2], :w1 = cos_pitch_angle(:u1, :B))
        @transform!(:Δw = :w1 .- :w0,)
        transform!([:α, :β] .=> ByRow(Int ∘ round ∘ rad2deg), renamecols=false)
        select!(Not(:wϕ0))
    end
end

function split_results(df)
    ldf = @subset(df, :sign .== 1)
    rdf = @subset(df, :sign .== -1)
    return ldf, rdf
end

subset_v(df) = @subset(df, :v .<= 128)
subset_β(df; βs=βs[2:2:end]) = @rsubset(df, :β ∈ βs)

#%%
subset_leave(df; tmax=tspan[2]) = @subset(df, :t1 .!= tmax)

rename_value_deg(v, s; kwargs...) = rename_value(v, s; kwargs...) * "°"
rename_value(v, s::Symbol; kwargs...) = "$(s) = $(round(v, kwargs...))"

rename_value_pair(v, s; r = rename_value, kwargs...) = v => r(v, s; kwargs...)
rename_value_pair_deg(v, s; r = rename_value_deg, kwargs...) = v => r(v, s; kwargs...)
vals(df, s) = unique(df[!, s]) |> sort


sign_map = :sign => renamer([1 => "Left-handed", -1 => "Right-handed"])
Δw_map = :Δw

function get_map(df, s; rename_func = rename_value_pair_deg)

    s => renamer(rename_func.(vals(df, s), s))
end

#%%
function pa_pair_plot(layer; axis=(;), kwargs...)
    v = visual(alpha=0.1, markersize=3; legend=(; markersize=12))
    l = layer * mapping(:w0, :w1) * v
    draw(l; axis, kwargs...)
end

w_axis = (; limits=((0, 1), (-1, 1)))

pa_pair_plot(df::AbstractDataFrame) = pa_pair_plot(data(df))

# Function to plot one figure per column
function pa_pair_plot_by_col!(df::AbstractDataFrame, s, axs, layer; scale=scales(), axis=w_axis, kwargs...)
    vs = (sort ∘ unique)(df[!, s])  # Get unique velocity values
    for (fg, v) in zip(axs, vs)
        df_s = @subset(df, $s .== v)
        plt = data(df_s) * mapping(:w0, :w1) * layer
        grids = draw!(fg, plt, scale; axis, kwargs...)

        r =  s == :v ? rename_value : rename_value_deg

        Label(fg[0, :], r(v, s))
        colorbar_pos = fg[1:size(grids, 1), size(grids, 2)+1]
        colorbar!(colorbar_pos, grids; scale=log10)
    end
end

pa_pair_plot_by_v!(df, axs) = pa_pair_plot_by_col!(df, :v, axs, mapping(col=α_map, row=β_map) * density_layer; scale)
pa_pair_plot_by_β!(df, axs) = pa_pair_plot_by_col!(df, :β, axs, mapping(col=α_map, row=v_map) * density_layer; scale)


"""
Plot the distribution of pitch-angle cosine variation
"""
function pa_diff_plot!(fig, layer; kwargs...)
    v = AlgebraOfGraphics.density() * visual(Lines)
    plt = layer * mapping(:Δw) * v
    axis = (; xlabel="Δw", ylabel="f (Δw)", limits=((-2, 1), (0, 4)))
    draw!(fig, plt; axis, kwargs...)
end

function pa_diff_plot(layer; kwargs...)
    v = AlgebraOfGraphics.density() * visual(Lines)
    plt = layer * mapping(:Δw) * v

    axis = (; xlabel="Δw", ylabel="f (Δw)", limits=((-2, 1), (0, 4)))
    draw(plt; axis, kwargs...)
end

pa_diff_plot!(fig, df::AbstractDataFrame; kwargs...) = pa_diff_plot!(fig, data(df) * mapping(col=α_map, row=β_map, color=v_map); kwargs...)
pa_diff_plot(df::AbstractDataFrame; kwargs...) = pa_diff_plot(data(df) * mapping(col=α_map, row=β_map, color=v_map); kwargs...)

begin
    density_layer = AlgebraOfGraphics.density(npoints=32)
    scale = scales(Color=(;
        colormap=:deep,
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
    plt = l * mapping(:w0, :ϕ0, color=color)
    draw(plt, scale; kwargs...)
end