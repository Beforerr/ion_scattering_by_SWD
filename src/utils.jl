using StaticArrays
using DrWatson

function extract_info(sol)
    u0 = SVector(sol.prob.u0...)
    u1 = SVector(sol.u[end]...)
    t1 = sol.t[end]
    return @dict u0 u1 t1
end