function m2(jumps)
    jumps_s = stack(jumps)
    diffs = jumps_s .- jumps_s[1, :]'
    m1 = mean(diffs, dims=2)
    m2 = mean(diffs .^ 2, dims=2) - m1 .^ 2
    vec(m2)
end