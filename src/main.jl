using DrWatson
using DataFrames
using DataFramesMeta
using CurrentSheetTestParticle

include("io.jl")
include("plot.jl")

params = [:θ, :β, :v, :sign]

subset_v(df) = @subset(df, :v .<= 128)
subset_v(v::Number) = df -> @subset(df, :v .<= v)
subset_β(itr) = df -> @rsubset(df, :β ∈ itr)

#%%
vals(df, s) = unique(df[!, s]) |> sort
rename_sym(s, v) = "$(s) = $(round(Int, v))"
rename_sym(s) = v -> rename_sym(s, v)
rename_sym_deg(s) = v -> (rename_sym(s, v) * "°")
get_map(s; renamer = rename_sym) = s => renamer(s)

θ_map = get_map(:θ; renamer = rename_sym_deg)
β_map = get_map(:β; renamer = rename_sym_deg)
v_map = get_map(:v)
sign_map = :sign => renamer([1 => "Left-handed", -1 => "Right-handed"])
w0 = :w0 => "w₀"
w1 = :w1 => "w₁"
α0 = :α0 => "α₀"
α1 = :α1 => "α₁"
Δα = :Δα
ϕ0 = :ϕ0 => "ϕ₀"
xyα = (α0, α1)
xyw = (w0, w1)


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

# Create Two-dimensional maps of the final value w1 are plotted as functions of the initial pitch-angle cosine w0 and gyrophase φ0, for α = 45◦ and β = 15◦ (upper panel), β = 45◦ (middle panel), β = 90◦ (lower panel).
function w1_map_plot(l; color = w1, scale= scales(;), kwargs...)
    plt = l * mapping(w0, ϕ0, color=color)
    draw(plt, scale; kwargs...)
end