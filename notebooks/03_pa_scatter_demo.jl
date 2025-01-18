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
tms = d["tms"]
ergs_approx = ["~10 eV", "~100 eV", "~5 keV", "~100 keV", "~1 MeV"]

@memoize matmul(tm, n::Int) = n == 0 ? I : tm * matmul(tm, n - 1)

"""
First moment would approximate zero
"""
function mix_rate(tm;)
    numBins = size(tm, 1)
    shiftedStates = [i - j for i in 1:numBins, j in 1:numBins] / numBins
    shiftedStatesS2 = shiftedStates .^ 2
    p = fill(1 / numBins, numBins)

    return n -> begin
        tmn = matmul(tm, n)
        m1 = p ⋅ diag(tmn * shiftedStates) # first moment
        m2 = p ⋅ diag(tmn * shiftedStatesS2) - m1^2
        return m2
    end
end

r = 0:200
res = map(tms) do tm
    mix_rate(tm).(r)
end

fit_label(f) = L"D_{nn} = %$(round(1 / f.b, digits=2))"
# Plot the mixing rate
function plot_mixing_rate!(res)
    options = Options(fine=N, nbest=8, besttol=1e-5)
    map(res) do re
        N = length(re)
        sca = scatter!(1:N, re)
        fit = fitexp(1:N, re; options)
        lin = lines!(fit.x, fit.y)
        [sca, lin, fit]
    end
end

let labels = ergs_approx, ax = (xlabel="n", ylabel="Mixing rate", xscale=log10, yscale=log10)
    f = Figure()
    ax = Axis(f[1, 1]; ax...)
    contents = plot_mixing_rate!(res[2:end])
    labels = LaTeXString.(ergs_approx[2:end] .* "(" .* fit_label.(getindex.(contents, 3)) .* ")")
    contents = getindex.(contents, Ref([1, 2]))
    axislegend(ax, contents, labels, "Energy"; position=:rb)
    xlims!(1.5, 200)
    ylims!(3e-2, 3e-1)
    easy_save("mixing_rate")
end

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
