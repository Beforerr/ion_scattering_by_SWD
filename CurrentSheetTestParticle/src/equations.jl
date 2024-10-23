using TestParticle: FTLError
using TestParticle: c

const c2 = c^2

"""
Normalized ODE equations for charged particle moving in static magnetic field with in-place form.

Optimized version of `trace_normalized!` for the case of static magnetic field.
"""
function trace_normalized_B!(du, u, p, t)
    _, _, B = p

    vx, vy, vz = @view u[4:6]
    Bx, By, Bz = B(u, t)

    du[1], du[2], du[3] = vx, vy, vz
    du[4] = vy*Bz - vz*By
    du[5] = vz*Bx - vx*Bz
    du[6] = vx*By - vy*Bx
end

"""
Normalized ODE equations for relativistic charged particle moving in static magnetic field with in-place form.
"""
function trace_relativistic_normalized_B!(du, u, p, t)
    q2m, B, v0 = p

    vx, vy, vz = @view u[4:6]
    Bx, By, Bz = B(u, t)

    u2 = (vx^2 + vy^2 + vz^2) * v0^2
    if u2 ≥ c2
        throw(DomainError(u2, FTLError))
    end
    γInv = √(1.0 - u2/c2)
 
    du[1], du[2], du[3] = vx, vy, vz
    du[4] = q2m * γInv *(vy*Bz - vz*By)
    du[5] = q2m * γInv *(vz*Bx - vx*Bz)
    du[6] = q2m * γInv *(vx*By - vy*Bx)
end