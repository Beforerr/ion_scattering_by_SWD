using OrdinaryDiffEq
using LinearAlgebra
using StaticArrays

include("fieldline.jl")
include("plot.jl")
include("utils.jl")
include("pa.jl")

function RD_B_field(r, α, β, B=1)
    z = r[3]
    ψ = β * tanh(z)
    Bz = B * cos(α)
    Bx = B * sin(α) * cos(ψ)
    By = B * sin(α) * sin(ψ)
    return SVector(Bx, By, Bz)
end