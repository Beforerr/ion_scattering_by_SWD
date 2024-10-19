# output the results as transition matrix
using DrWatson
@quickactivate
using FHist
using Comonicon: @main, @cast
include("../src/io.jl")

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
    @unpack dir, bins = d
    w_range = range(-1, 1, length=bins+1)
    df = get_result_dfs(; dir)
    df.tm = map(eachrow(df)) do r
        transition_matrix_w(get_result(r), w_range)
    end
    select!(df, Not(:results))
    @dict df w_range
end

"""
Calculate the transition matrix for each simulation result and save it in a file

# Options

- `-b, --bins`: Number of bins for the transition matrix
- `--dir`: directory to read the results from
"""
@main function main(; bins::Int = 45, dir::String = "simulations", path=datadir(), prefix="tm")
    d = @dict bins dir
    produce_or_load(produce, d, path; prefix)
end