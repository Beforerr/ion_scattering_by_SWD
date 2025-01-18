const DEFAULT_θ = 45
const DEFAULT_β = 90
const DEFAULT_SIGN = 1


"""
    RotationDiscontinuity

Rotating magnetic field. 

# Arguments
- `θ` : azimuthal angle, angle between Bn and B0
- `β` : half of rotation angle
- `sign` : sign of rotation, 1 for left-handed, -1 for right-handed

# Notes
φ = β * tanh(z) is the polar angle
"""
@kwdef struct RotationDiscontinuity
    B = 1
    θ = DEFAULT_θ
    β = DEFAULT_β
    sign::Integer = DEFAULT_SIGN
end

"""
# Keyword Arguments
- `dir` : direction where the field depends on, 1, 2, or 3
"""
function B(r, conf::RotationDiscontinuity; dir=3)
    @unpack θ, β, sign, B = conf
    z = r[dir]
    φ = β * tanh(z)
    Bz = B * cosd(θ)
    Bx = B * sind(θ) * sind(φ)
    By = sign * B * sind(θ) * cosd(φ)
    return SVector(Bx, By, Bz)
end

B(conf::RotationDiscontinuity; kwargs...) = r -> B(r, conf; kwargs...)

function TD_B_field(r; dir=3, B=1, By=0, θ=DEFAULT_θ, β=DEFAULT_β, kw...)
    z = r[dir]
    φ = β * tanh(z)
    Bz = B * cosd(θ)
    Bx = B * sind(θ) * sind(φ)
    return SVector(Bx, By, Bz)
end

TD_B_field(; kw...) = r -> TD_B_field(r; kw...)