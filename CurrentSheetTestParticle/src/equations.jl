using TestParticle: FTLError
using TestParticle: c
using LinearAlgebra: ×

const c2 = c^2

"""
Reciprocal of Lorentz factor for relativistic particle.
"""
function gamma_reciprocal(vx, vy, vz; v0=1)
    u2 = (vx^2 + vy^2 + vz^2) * v0^2
    u2 ≥ c2 && throw(DomainError(u2, FTLError))
    return √(1 - u2 / c2)
end

"""
Normalized ODE equations for charged particle moving in static magnetic field with in-place form.

Optimized version of `trace_normalized!` for the case of static magnetic field.
"""
function trace_normalized_B!(du, u, p, t)
    _, _, B = p

    v = @views SVector{3}(u[4:6])
    b = SVector{3}(B(u, t))

    du[1:3] = v
    du[4:6] = v × b
end

"""
Out-of-place version of `trace_normalized_B!`, with a little bit better performance.
"""
function trace_normalized_B(u, p, t)
    _, _, B = p

    v = @views SVector{3}(u[4:6])
    b = SVector{3}(B(u, t))

    dv = v × b
    return SVector{6}(v..., dv...)
end

function trace_normalized_B_1D!(du, u, p, t; dir=1)
    _, _, B = p
    v = @views SVector{3}(u[4:6])
    b = SVector{3}(B(u, t))
    du[1] = v[dir]
    du[2:4] = v × b
end