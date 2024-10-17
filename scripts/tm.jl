# output the results as transition matrix
using DrWatson
@quickactivate
using Revise
using FHist
using ArgParse
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
    @unpack dir, bins = d
    w_range = range(-1, 1, length=bins+1)
    df = get_result_dfs(; dir)
    df.tm = map(eachrow(df)) do r
        transition_matrix_w(get_result(r), w_range)
    end
    select!(df, Not(:results, :diffeq))
    @dict df w_range
end

function main(; path=datadir(), prefix="tm")
    parsed_args = parse_commandline()
    produce_or_load(produce, parsed_args, path; prefix)
end

function parse_commandline(; dir = "test_alg", bins = 45)
    s = ArgParseSettings()
    @add_arg_table! s begin
        "--dir", "-d"
        help = "directory to read the results from"
        default = dir
        "--bins", "-b"
        help = "Number of bins for the transition matrix"
        arg_type = Int
        default = bins
    end
    return parse_args(s)
end

main()