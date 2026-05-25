module Zinq

using FFTW, Dates, LinearAlgebra, Printf, Tullio

export julia_main

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

function get_prop_k(m::Float64, k::NTuple{N, AbstractArray{Float64}}, dt::Float64) where N
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

function print_iter(i::Int, pe::Float64, ke::Float64, pos::NTuple{N, Float64}, mom::NTuple{N, Float64}, pop::NTuple{S, Float64}, norm::Float64, elapsed::Float64) where {N, S}
    fmt_pos = join([@sprintf("%10.4f", x) for x in pos], " ")
    fmt_mom = join([@sprintf("%10.4f", x) for x in mom], " ")
    fmt_pop = join([@sprintf("%10.4f", x) for x in pop], " ")

    values = (i - 1, ke, pe, ke + pe, fmt_pos, fmt_mom, fmt_pop, norm, format_duration(elapsed))

    @printf("%7d %12.6f %12.6f %12.6f [%s] [%s] [%s] %9.4f %s\n", values...)
end

function print_header(ndim::Int, nstate::Int)
    dim_w, state_w = 11 * ndim + 1, 11 * nstate + 1

    labels = ("ITER", "KIN (Eh)", "POT (Eh)", "TOT (Eh)", dim_w, "POS (a0)", dim_w, "MOM (hb/a0)", state_w, "POPULATION", "NORM", "TIME")

    @printf("%7s %12s %12s %12s %*s %*s %*s %9s %-s\n", labels...)
end

function format_duration(seconds::Float64)
    return Dates.format(Time(0) + Millisecond(round(Int, seconds * 1000)), "HH:MM:SS.sss")
end

function run_qd()
    ic = InitialConditions((-10.0,), (15.0,), (2.0,), 1, true)
    grid = Grid(((-24.0, 32.0),), 2048)
    m = 2000.0
    iters = 3000
    dt = 1.0
    log_interval = 500

    r = gen_grid_r(grid)
    k = gen_grid_k(grid)

    W = gen_wfn(ic, r, grid.npoint, 2)
    V = tully_1(r..., 0.01)
    A, U = get_pot_eigen(V)
    K = get_prop_k(m, k, dt)
    R = get_prop_r(A, U, dt / 2)

    normalize!(W, grid)

    print_header(length(r), size(W, ndims(W)))

    for i in 1:iters + 1
        elapsed = @elapsed if i > 1
            propagate!(W, R, K)
        end

        if (i - 1) % log_interval == 0
            W_k = to_kspace(W)

            norm = calc_norm(W, grid)
            pe = calc_pe(W, V, grid)
            ke = calc_ke(W_k, grid, m, k)
            pop = calc_pop(W, grid)
            pos = calc_pos(W, grid, r)
            mom = calc_mom(to_kspace(W), grid, k)

            print_iter(i, pe, ke, pos, mom, pop, norm, elapsed)
        end
    end
end

function julia_main()
    run_qd()
end

end # module Zinq
