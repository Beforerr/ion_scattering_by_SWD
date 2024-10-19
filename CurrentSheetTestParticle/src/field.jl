const DEFAULT_θ = 45
const DEFAULT_β = 90
const DEFAULT_SIGN = 1

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
function RD_B_field(r; B=1, θ=DEFAULT_θ, β=DEFAULT_β, sign=DEFAULT_SIGN)
    z = r[3]
    φ = β * tanh(z)
    Bz = B * cosd(θ)
    Bx = B * sind(θ) * sind(φ)
    By = sign * B * sind(θ) * cosd(φ)
    return SVector(Bx, By, Bz)
end

RD_B_field(; kwargs...) = r -> RD_B_field(r; kwargs...)