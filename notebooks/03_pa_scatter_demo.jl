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
    TM_OBS_FILE = "tm_obs.jld2"
end

# Load the data
d = wload(datadir(TM_OBS_FILE))
df = d["df"]
vPs = d["vPs"]
tms = d["tms"]
ergs_approx = ["~10 eV", "~100 eV", "~5 keV", "~100 keV"]

@memoize matmul(tm, n::Int) =  n == 0 ? I : tm * matmul(tm, n-1)

"""
First moment would approximate zero
"""
function mix_rate(tm)
    numBins = size(tm, 1)
    shiftedStates = [i-j for i in 1:numBins, j in 1:numBins]
    shiftedStatesS2 = shiftedStates .^2
    p = fill(1/numBins, numBins)

    return n -> begin
        tmn = matmul(tm, n) 
        m1 = p ⋅ diag(tmn * shiftedStates) # first moment
        m2 = p ⋅ diag(tmn * shiftedStatesS2) # second moment
        m2 - m1^2
    end
end

r = 0:100
res = map(tms) do tm
    mix_rate(tm).(r)
end

# Plot the mixing rate
let labels = ergs_approx, ax = (xlabel="n", ylabel="Mixing rate", xscale=log10, yscale=log10)
    f = Figure()
    ax = Axis(f[1,1]; ax...)
    contents = stairs!.(res)
    Legend(f[1, 2], contents, labels, "Energy")
    easy_save("mixing_rate")
end

α

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
v0s = vs
w0s = sample(ws, length(vs))

# Iterate the jump n times and plot the 
n = 100
w_hist = pa_jump.(w0s, v0s; n)

let labels = ergs_approx
    f = Figure()
    ax = Axis(f[1,1], xlabel="n", ylabel="Cos(α)")
    contents = stairs!.(w_hist)
    Legend(f[1, 2], contents, labels, "Energy")
    easy_save("pa_jump_history")
end
