function m2(jumps)
    jumps_s = stack(jumps)
    diffs = jumps_s .- jumps_s[1, :]'
    m1 = mean(diffs, dims=2)
    m2 = mean(diffs .^ 2, dims=2) - m1 .^ 2
    vec(m2)
end


@memoize matmul(tm, n::Int) = n == 0 ? I : tm * matmul(tm, n - 1)

"""
Mixing rates using transition matrix
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