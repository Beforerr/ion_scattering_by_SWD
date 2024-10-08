"""
Cos pitch angle
"""
function cos_pitch_angle(u, B)
    v = @view u[4:6]
    return B ⋅ v / norm(v) / norm(B)
end
    
function cos_pitch_angle(sol::ODESolution, i)
    B_field = r -> RD_B_field(r, α, β)
    B_func = sol.prob.p[3]
    u = sol.u[i]
    B = B_func(u)
    return cos_pitch_angle(u, B)
end