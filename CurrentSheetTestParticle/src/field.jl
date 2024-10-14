"""
    RD_B_field([r]; ...)

Rotating magnetic field. 

# Keyword Arguments
- `θ` : azimuthal angle, angle between Bn and B0
- `β` : half of rotation angle
- `sign` : sign of rotation, 1 for left-handed, -1 for right-handed

# Notes
φ = β * tanh(z) is the polar angle
"""
function RD_B_field(r; B=1, θ=π/2, β=π/2, sign=1)
    z = r[3]
    φ = β * tanh(z)
    Bz = B * cos(θ)
    Bx = B * sin(θ) * sin(φ)
    By = sign * B * sin(θ) * cos(φ)
    return SVector(Bx, By, Bz)
end

RD_B_field(; kwargs...) = r -> RD_B_field(r; kwargs...)