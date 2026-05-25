module Potential

export POTENTIALS

struct Potential
    fn::Function
    nstate::Int
end

function tully_1(r::AbstractArray{Float64}; A::Float64 = 0.01, B::Float64 = 1.6, C::Float64 = 0.005, D::Float64 = 1.0)
    V = Array{Float64}(undef, 2, 2, size(r)...)

    @. V[1, 1, :] = sign(r) * A * (1 - exp(-B * abs(r)))
    @. V[1, 2, :] = C * exp(-D * r^2)
    @. V[2, 1, :] = V[1, 2, :]
    @. V[2, 2, :] = -V[1, 1, :]

    return V
end

function harmonic(r::Vararg{AbstractArray{Float64}, N}; k::NTuple{N, Float64}=ntuple(_ -> 1.0, N)) where N
    return reshape(sum(i -> 0.5 * k[i] .* r[i].^2, 1:N), 1, 1, size(r[1])...)
end

const POTENTIALS = Dict{String, Potential}(
    "harmonic" => Potential(harmonic, 1),
    "tully_1" => Potential(tully_1, 2),
)

end # module Potential
