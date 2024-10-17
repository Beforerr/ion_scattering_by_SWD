# output the results as transition matrix
using DrWatson
@quickactivate
using Revise
using FHist
using CurrentSheetTestParticle
using JLD2
includet("../src/io.jl")

"""
right stochastic matrix, with each row summing to 1., m〚i,j〛 = Probability[x[k+1]=j|x[k]=i]
"""
function transition_matrix(h::Hist2D)
    m = bincounts(h)
    m ./= sum(m, dims=2)
end

function transition_matrix_w(df, w_range)
    h = Hist2D((df.w0, df.w1); binedges=(w_range, w_range))
    transition_matrix(h)
end

function produce(d::Dict)
    @unpack dir, w_range = d
    df = collect_results(datadir(dir))
    df.tm = map(eachrow(df)) do r
        transition_matrix_w(subset_leave(get_result(r)), w_range)
    end
    select!(df, Not(:results, :diffeq))
    @dict df w_range
end


path = datadir()
params = dict_list(Dict(
    :dir => ["test_alg"],
    :w_range => range(-1, 1, length=46)
))

param_dicts = map(params) do d
    produce_or_load(produce, d, path; prefix = "tm")
end