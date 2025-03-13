begin
    using DrWatson
    @quickactivate
    using StatsBase
    using StatsBase: sample
    using DataFrames
    using Unitful
    using Beforerr
    using DimensionalData
    using LinearAlgebra
    using Memoization
    using CairoMakie
    using LaTeXStrings
    using EasyFit
    TM_OBS_FILE = "tm_obs.jld2"

    update_theme!(Beforerr.theme_pub())
end

# Load the data
d = wload(datadir(TM_OBS_FILE))
df = d["df"]
vPs = d["vPs"]
ergs_approx = ["~10 eV", "~100 eV", "~5 keV", "~100 keV", "~1 MeV"]

"""
Pitch angle jump

given the w0, sample w1 with probability proportional to the value
"""
function pa_jump(w0, v)
    df_s = @subset(df, :w0 .== w0, :v .== v)
    w1s, probs = df_s.w1, df_s.value
    return sample(w1s, Weights(probs))
end

function pa_jump(w0, args...; n)
    w = zeros(n)
    w[1] = w0
    for i in 2:n
        w[i] = pa_jump(w[i-1], args...)
    end
    return w
end

ws = unique(df.w0) |> sort
vs = unique(df.v) |> sort

begin
    v0s = vs
    w0s = sample(ws, length(vs))
    n = 100
    w_hist = pa_jump.(w0s, v0s; n)
end


# plot the history
let labels = ergs_approx
    f = Figure()
    ax = Axis(f[1, 1], xlabel="n", ylabel="Cos(α)")
    contents = stairs!.(w_hist)
    Legend(f[1, 2], contents, labels, "Energy")
    easy_save("pa_jump_history")
end


let labels = ergs_approx, idx = 4:5
    f = Figure()
    ax = Axis(f[1, 1], xlabel="n", ylabel="Cos(α)")
    contents = stairs!.(w_hist[idx])
    Legend(f[1, 2], contents, labels[idx], "Energy")
    easy_save("pa_jump_history_high")
end
