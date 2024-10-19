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
function RD_B_field(r; B=1, θ=45, β=90, sign=1)
    z = r[3]
    φ = β * tanh(z)
    Bz = B * cosd(θ)
    Bx = B * sind(θ) * sind(φ)
    By = sign * B * sind(θ) * cosd(φ)
    return SVector(Bx, By, Bz)
end

RD_B_field(; kwargs...) = r -> RD_B_field(r; kwargs...)