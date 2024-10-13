"""
Rotating magnetic field. 

# Notes
Left-handed rotating with `sign=1`.
Right-handed rotating with `sign=-1`
"""
function RD_B_field(r, α, β; B=1, sign=1)
    z = r[3]
    ψ = β * tanh(z)
    Bz = B * cos(α)
    Bx = B * sin(α) * sin(ψ)
    By = sign * B * sin(α) * cos(ψ)
    return SVector(Bx, By, Bz)
end

RD_B_field(; B=1, α, β, sign) = r -> RD_B_field(r, α, β; B, sign)