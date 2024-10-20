begin
	using DrWatson
    @quickactivate
    using StatsBase
	using StatsBase: sample
	using DataFrames
	using Unitful
    using Beforerr
    
    TM_OBS_FILE = "tm_obs.jld2"
end

# Load the data
df = wload(datadir(TM_OBS_FILE))["df"]
ws = unique(df.w0) |> sort
vs = unique(df.v) |> sort

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

v0s = vs
w0s = sample(ws, length(vs))

# Iterate the jump n times and plot the 
n = 100
w_hist = pa_jump.(w0s, v0s; n)

# Plot the jump history using Stairs
using CairoMakie

ergs_approx = ["~10 eV", "~100 eV", "~5 keV", "~100 keV"]

let labels = ergs_approx
    f = Figure()
    ax = Axis(f[1,1], xlabel="n", ylabel="Cos(α)")
    contents = stairs!.(w_hist)
    Legend(f[1, 2], contents, labels, "Energy")
    easy_save("pa_jump_history")
end
