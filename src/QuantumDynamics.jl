module QuantumDynamics

include("Potential.jl")

using .Potential, FFTW, Dates, LinearAlgebra, Printf, TimerOutputs, Tullio

export run_qd

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

struct SimulationContext{N}
    W::AbstractArray{ComplexF64}
    R::AbstractArray{ComplexF64}
    K::AbstractArray{ComplexF64}
    V::AbstractArray{Float64}
    A::AbstractArray{Float64}
    U::AbstractArray{Float64}
    r::NTuple{N, AbstractArray{Float64}}
    k::NTuple{N, AbstractArray{Float64}}
    m::Float64
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

function gen_wfn(ic::InitialConditions{N}, r::NTuple{N, AbstractArray{Float64}}, npoint::Int, nstate::Int) where N
    W = zeros(ComplexF64, (ntuple(_ -> npoint, N)..., nstate))

    components = ntuple(N) do i
        @. exp(-0.5ic.gamma[i] * (r[i] - ic.pos[i])^2 + im * ic.mom[i] * (r[i] - ic.pos[i]))
    end

    broadcast!(*, selectdim(W, N + 1, ic.state + 1), components...)

    return W
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

function get_prop_r(A::AbstractArray{Float64}, U::AbstractArray{Float64}, dt::ComplexF64)
    R = similar(U, ComplexF64)

    for l in CartesianIndices(size(A)[2:end])
        @views R[:, :, l] = U[:, :, l] * Diagonal(@. exp(-im * A[:, l] * dt)) * U[:, :, l]'
    end

    return R
end

function get_prop_k(m::Float64, k::NTuple{N, AbstractArray{Float64}}, dt::ComplexF64) where N
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

function print_iter(i::Int, obs::Observables{N}, elapsed::UInt64) where N
    fmt_pos = join([@sprintf("%10.4f", x) for x in obs.pos], " ")
    fmt_mom = join([@sprintf("%10.4f", x) for x in obs.mom], " ")
    fmt_pop = join([@sprintf("%10.4f", x) for x in obs.pop], " ")

    values = (i - 1, obs.ke, obs.pe, obs.ke + obs.pe, fmt_pos, fmt_mom, fmt_pop, obs.norm, format_duration(elapsed))

    @printf("%7d %12.6f %12.6f %12.6f [%s] [%s] [%s] %9.4f %s\n", values...)
end

function print_header(state::Int, ndim::Int, nstate::Int, itp::Bool)
    if itp > 0
        print("\nSTATE $state ITP")
    end

    dim_w, state_w = 11 * ndim + 1, 11 * nstate + 1

    labels = ("ITER", "KIN (Eh)", "POT (Eh)", "TOT (Eh)", dim_w, "POS (a0)", dim_w, "MOM (hb/a0)", state_w, "POPULATION", "NORM", "TIME")

    @printf("\n%7s %12s %12s %12s %*s %*s %*s %9s %-s\n", labels...)
end

function format_duration(nanos::UInt64)
    return Dates.format(Time(0) + Nanosecond(nanos), "HH:MM:SS.sss")
end

function to_adia(W::AbstractArray{ComplexF64}, U::AbstractArray{Float64})
    W_a, W_f = similar(flat_wfn(W)), flat_wfn(W)

    @tullio W_a[i, I] := U[I, J, i]' * W_f[i, J]

    return reshape(W_a, size(W))
end

function calc_observables(ctx::SimulationContext{N}, grid::Grid{N}, adia::Bool) where N
    W_k = to_kspace(ctx.W)

    norm = calc_norm(ctx.W, grid)
    pe = calc_pe(ctx.W, ctx.V, grid)
    ke = calc_ke(W_k, grid, ctx.m, ctx.k)
    pop = calc_pop(adia ? to_adia(ctx.W, ctx.U) : ctx.W, grid)
    pos = calc_pos(ctx.W, grid, ctx.r)
    mom = calc_mom(W_k, grid, ctx.k)

    return Observables(pos, mom, pop, norm, pe, ke)
end

function overlap(W1::AbstractArray{ComplexF64}, W2::AbstractArray{ComplexF64}, grid::Grid{N}) where N
    return sum(conj(W1[l]) * W2[l] for l in CartesianIndices(W1)) * get_dr(grid)
end

function project_out!(W1::AbstractArray{ComplexF64}, W2::AbstractArray{ComplexF64}, grid::Grid{N}) where N
    W1 .-= overlap(W2, W1, grid) .* W2
end

function log_final_pop(W::AbstractArray{ComplexF64}, grid::Grid{N}) where N
    pop = calc_pop(W, grid)

    for (i, p) in enumerate(pop)
        @printf("%sFINAL POPULATION OF STATE %02d: %.6f\n", i == 1 ? "\n" : "", i, p)
    end
end

function log_final_te(tes::Vector{Float64})
    for (i, te) in enumerate(tes)
        @printf("%sENERGY OF STATE %02d: %8.6f Eh\n", i == 1 ? "\n" : "", i, te)
    end
end

function project_and_normalize!(W::AbstractArray{ComplexF64}, wfns::Vector{AbstractArray{ComplexF64}}, grid::Grid{N}) where N
    for W2 in wfns
        project_out!(W, W2, grid)
    end

    normalize!(W, grid)
end

function run_qd()
    # ic = InitialConditions((-10.0,), (15.0,), (2.0,), 1, true)
    # grid = Grid(((-24.0, 32.0),), 4096)
    # m = 2000.0
    # iters = 3000
    # dt = 1.0
    # log_interval = 500
    # adia = true
    # itp = false
    # pot = POTENTIALS["tully_1"]

    ic = InitialConditions((1.0,), (0.0,), (2.0,), 0, true)
    grid = Grid(((-8.0, 8.0),), 256)
    m = 1.0
    iters = 1000
    dt = 0.01
    log_interval = 200
    adia = false
    itp = 5
    pot = POTENTIALS["harmonic"]

    @timeit "INITIALIZATION" begin
        r = gen_grid_r(grid)
        k = gen_grid_k(grid)

        V = pot.fn(r...); A, U = get_pot_eigen(V)

        R = get_prop_r(A, U, 0.5 * (itp > 0 ? -im * dt : dt + 0im))
        K = get_prop_k(m, k, 1.0 * (itp > 0 ? -im * dt : dt + 0im))
    end

    opt_wfn, opt_wfn_te = AbstractArray{ComplexF64}[], Float64[]

    for i in 1:(itp > 0 ? itp : 1)

        W = gen_wfn(ic, r, grid.npoint, pot.nstate); normalize!(W, grid)

        ctx = SimulationContext(W, R, K, V, A, U, r, k, m)

        print_header(i, length(r), pot.nstate, itp > 0)

        start_time = time_ns()

        for j in 1:iters + 1
            if j > 1
                @timeit "PROPAGATION" propagate!(W, R, K)
            end

            @timeit "PROJECTION" if itp > 0
                project_and_normalize!(W, opt_wfn, grid)
            end

            if (j - 1) % log_interval == 0
                @timeit "OBSERVABLES" obs = calc_observables(ctx, grid, adia)

                print_iter(j, obs, time_ns() - start_time)

                if j == iters + 1 && itp > 1
                    push!(opt_wfn_te, obs.pe + obs.ke)
                end

                start_time = time_ns()
            end
        end

        if pot.nstate > 1 && itp == 0
            log_final_pop(W, grid)
        end

        push!(opt_wfn, W)
    end

    if itp > 1
        log_final_te(opt_wfn_te)
    end
end

end # module QuantumDynamics
