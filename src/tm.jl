using FHist

const DEFAULT_W_BINEDGE = range(-1, 1, length=46)
const DEFAULT_W_BINEDGES = (DEFAULT_W_BINEDGE, DEFAULT_W_BINEDGE)

"""
right stochastic matrix, with each row summing to 1., m〚i,j〛 = Probability[x[k+1]=j|x[k]=i]
"""
function transition_matrix(h::Hist2D)
    m = bincounts(h)
    m ./= sum(m, dims=2)
end

function transition_matrix_w(df; binedges=DEFAULT_W_BINEDGES, kw...)
    h = Hist2D((df.w0, df.w1); binedges, kw...)
    transition_matrix(h)
end