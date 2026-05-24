module Zinq

using FFTW, LinearAlgebra, Tullio

export Grid, InitialConditions, gen_grid_r, gen_grid_k, gen_wfn, tully_1, normalize!, calc_norm, calc_pe, calc_pop, calc_pos, to_kspace, calc_ke, calc_mom, get_pot_eigen, get_prop_r, get_prop_k, propagate!

struct Grid{N}
    limits::NTuple{N, Tuple{Float64, Float64}}
    npoint::Int
end

struct InitialConditions{N}
    pos::NTuple{N, Float64}
    mom::NTuple{N, Float64}
    gamma::NTuple{N, Float64}
    state::Int
    adia::Bool
end

struct Observables{N, S}
    pos::NTuple{N, Float64}
    mom::NTuple{N, Float64}
    pop::NTuple{S, Float64}
    norm::Float64
    pe::Float64
    ke::Float64
end

function gen_grid_r(grid::Grid{N}) where N
    limits, npoint = grid.limits, grid.npoint

    axes = ntuple(N) do i
        Vector(range(limits[i][1], step=(limits[i][2] - limits[i][1]) / npoint, length=npoint))
    end

    return ntuple(N) do i
        reshape(axes[i], ntuple(j -> j == i ? (:) : 1, N)...)
    end
end

function gen_grid_k(grid::Grid{N}) where N
    limits, npoint = grid.limits, grid.npoint

    axes = ntuple(N) do i
        Vector(fftfreq(npoint, 2pi * npoint / (limits[i][2] - limits[i][1])))
    end

    return ntuple(N) do i
        reshape(axes[i], ntuple(j -> j == i ? (:) : 1, N)...)
    end
end

function gen_wfn(ic::InitialConditions{N}, r::NTuple{N, AbstractArray{Float64}} , npoint::Int, nstate::Int) where N
    W = zeros(ComplexF64, (ntuple(_ -> npoint, N)..., nstate))

    components = ntuple(N) do i
        @. exp(-0.5ic.gamma[i] * (r[i] - ic.pos[i])^2 + im * ic.mom[i] * (r[i] - ic.pos[i]))
    end

    broadcast!(*, selectdim(W, N + 1, ic.state + 1), components...)

    return W
end

function tully_1(r::AbstractVector{Float64}, A::Float64 = 0.1, B::Float64 = 1.6, C::Float64 = 0.005, D::Float64 = 1.0)
    V = Array{Float64}(undef, 2, 2, size(r)...)

    @views begin
        @. V[1, 1, :] = sign(r) * A * (1 - exp(-B * abs(r)))
        @. V[1, 2, :] = C * exp(-D * r^2)
        @. V[2, 1, :] = V[1, 2, :]
        @. V[2, 2, :] = -V[1, 1, :]
    end

    return V
end

function get_dr(grid::Grid{N}) where N
    return prod(map(limit -> (limit[2] - limit[1]) / grid.npoint, grid.limits))
end

function calc_norm(W::AbstractArray{ComplexF64}, grid::Grid{N}) where N
    return sqrt(sum(abs2, W) * get_dr(grid))
end

function flat_wfn(W::AbstractArray{ComplexF64})
    return reshape(W, :, size(W, ndims(W)))
end

function flat_pot(V::AbstractArray{Float64})
    return reshape(V, size(V, 1), size(V, 2), :)
end

function calc_pe(W::AbstractArray{ComplexF64}, V::AbstractArray{Float64}, grid::Grid{N}) where N
    W_f, V_f = flat_wfn(W), flat_pot(V)

    @tullio pe := conj(W_f[i, I]) * V_f[I, J, i] * W_f[i, J]

    return real(pe) * get_dr(grid)
end

function calc_pop(W::AbstractArray{ComplexF64}, grid::Grid{N}) where N
    W_f = flat_wfn(W)
    
    @tullio pop[I] := abs2(W_f[i, I])
    
    return Tuple(pop .* get_dr(grid))
end

function normalize!(W::AbstractArray{ComplexF64}, grid::Grid{N}) where N
    W ./= calc_norm(W, grid) 
end

function calc_pos(W::AbstractArray{ComplexF64}, grid::Grid{N}, r::NTuple{N, AbstractArray{Float64}}) where N
    return ntuple(N) do k
        sum(abs2(W[l]) * r[k][l[k]] for l in CartesianIndices(W)) * get_dr(grid)
    end
end

function to_kspace(W::AbstractArray{ComplexF64})
    return fft(W, 1:ndims(W) - 1)
end

function get_dk(grid::Grid{N}) where N
    return get_dr(grid) / grid.npoint^N
end

function calc_mom(W::AbstractArray{ComplexF64}, grid::Grid{N}, k::NTuple{N, AbstractArray{Float64}}) where N
    return ntuple(N) do i
        sum(abs2(W[l]) * k[i][l[i]] for l in CartesianIndices(W)) * get_dk(grid)
    end
end

function calc_ke(W::AbstractArray{ComplexF64}, grid::Grid{N}, m::Float64, k::NTuple{N, AbstractArray{Float64}}) where N
    return 0.5 * sum(abs2(W[l]) * sum(k[i][l[i]]^2 / m for i in 1:N) for l in CartesianIndices(W)) * get_dk(grid)
end

function get_pot_eigen(V::AbstractArray{Float64, P}) where P
    A, U = similar(V, (size(V, 1), size(V)[3:end]...)), similar(V)

    for l in CartesianIndices(size(V)[3:end])
        A[:, l], U[:, :, l] = eigen(Symmetric(view(V, :, :, l)))
    end
    
    return A, U
end

function get_prop_r(A::AbstractArray{Float64}, U::AbstractArray{Float64}, dt::Float64)
    R = similar(U, ComplexF64)

    for l in CartesianIndices(size(A)[2:end])
        @views R[:, :, l] = U[:, :, l] * Diagonal(@. exp(-im * A[:, l] * dt)) * U[:, :, l]'
    end

    return R
end

function get_prop_k(grid::Grid{N}, m::Float64, k::NTuple{N, AbstractArray{Float64}}, dt::Float64) where N
    return @. exp(-0.5im * $(reduce(.+, @. k[i]^2 / m for i in 1:N)) * dt)
end

function propagate_r!(W::AbstractArray{ComplexF64}, R::AbstractArray{ComplexF64})
    W_f = flat_wfn(W)

    @tullio W_t[i, I] := R[I, J, i] * W_f[i, J]

    W .= W_t
end

function propagate_k!(W::AbstractArray{ComplexF64}, K::AbstractArray{ComplexF64}) 
    fft!(W, 1:ndims(W) - 1)

    W .*= K

    ifft!(W, 1:ndims(W) - 1)
end

function propagate!(W::AbstractArray{ComplexF64}, R::AbstractArray{ComplexF64}, K::AbstractArray{ComplexF64})
    propagate_r!(W, R)
    propagate_k!(W, K)
    propagate_r!(W, R)
end

end # module Zinq
