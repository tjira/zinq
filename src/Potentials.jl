module Potentials

export POTENTIALS, Potential, get_nstate

struct Potential{S, F <: Function}
    fn::F; td::Bool
end

function harmonic(t::Float64, r::Vararg{AbstractArray{Float64}, N}; k::NTuple{N, Float64}=ntuple(_ -> 1.0, N)) where N
    V = @. $reduce(.+, 0.5 * k[i] * r[i]^2 for i in 1:N)

    return reshape(V, 1, 1, size(V)...)
end

function time_linear(t::Float64, r::Vararg{AbstractArray{Float64}, N}; coupling::Float64 = 2.0, slope::Float64 = 10.0) where N
    V = Array{Float64}(undef, 2, 2, size(r[1])...)

    @. V[1, 1, :] = slope * (t - slope)
    @. V[1, 2, :] = coupling
    @. V[2, 1, :] = coupling
    @. V[2, 2, :] = -V[1, 1, :]

    return V
end

function tully_1(t::Float64, r::Vararg{AbstractArray{Float64}, N}; A::Float64 = 0.01, B::Float64 = 1.6, C::Float64 = 0.005, D::Float64 = 1.0) where N
    V = Array{Float64}(undef, 2, 2, size(r[1])...)

    @. V[1, 1, :] = sign(r[1]) * A * (1 - exp(-B * abs(r[1])))
    @. V[1, 2, :] = C * exp(-D * r[1]^2)
    @. V[2, 1, :] = V[1, 2, :]
    @. V[2, 2, :] = -V[1, 1, :]

    return V
end

function get_nstate(::Potential{S, F}) where {S, F}
    return S
end

const POTENTIALS = Dict(
    "harmonic"    => Potential{1, typeof(harmonic   )}(harmonic,    false),
    "time_linear" => Potential{2, typeof(time_linear)}(time_linear, true ),
    "tully_1"     => Potential{2, typeof(tully_1    )}(tully_1,     false),
)

end # module Potentials
