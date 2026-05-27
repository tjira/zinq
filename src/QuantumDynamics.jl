module QuantumDynamics

include("Potentials.jl")

using .Potentials, FFTW, HDF5, TOML, Dates, LinearAlgebra, Printf, TimerOutputs, Tullio

export run_qd

Base.@kwdef struct Absorber{N}
    bounds::NTuple{N, Tuple{Float64, Float64}}; exponent::Float64
end

Base.@kwdef struct Grid{N}
    bounds::NTuple{N, Tuple{Float64, Float64}}; npoint::Int
end

Base.@kwdef struct InitialConditions{N}
    pos  ::NTuple{N, Float64}
    mom  ::NTuple{N, Float64}
    gamma::NTuple{N, Float64}

    state::Int; adia::Bool
end

Base.@kwdef struct Observables{N, S}
    pos::Union{NTuple{N, Float64}, Nothing}
    mom::Union{NTuple{N, Float64}, Nothing}
    pop::Union{NTuple{S, Float64}, Nothing}

    norm::Union{Float64, Nothing}
    pe  ::Union{Float64, Nothing}
    ke  ::Union{Float64, Nothing}
end

Base.@kwdef struct History{N, S}
    pos::Union{Vector{NTuple{N, Float64}}, Nothing} = nothing
    mom::Union{Vector{NTuple{N, Float64}}, Nothing} = nothing
    pop::Union{Vector{NTuple{S, Float64}}, Nothing} = nothing

    norm::Union{Vector{Float64}, Nothing} = nothing
    pe  ::Union{Vector{Float64}, Nothing} = nothing
    ke  ::Union{Vector{Float64}, Nothing} = nothing

    wfn::Union{Vector{AbstractArray{ComplexF64}}, Nothing} = nothing
end

Base.@kwdef struct Simulation
    itp         ::Int
    iters       ::Int
    log_interval::Int

    dt::Float64
    m ::Float64

    adia::Bool
end

Base.@kwdef struct Writer
    norm::Union{String, Nothing} = nothing
    pos ::Union{String, Nothing} = nothing
    mom ::Union{String, Nothing} = nothing
    pop ::Union{String, Nothing} = nothing
    pe  ::Union{String, Nothing} = nothing
    ke  ::Union{String, Nothing} = nothing
    te  ::Union{String, Nothing} = nothing
    wfn ::Union{String, Nothing} = nothing
end

Base.@kwdef struct SimulationContext{N}
    W::AbstractArray{ComplexF64}
    R::AbstractArray{ComplexF64}
    K::AbstractArray{ComplexF64}

    V::AbstractArray{Float64}
    A::AbstractArray{Float64}
    U::AbstractArray{Float64}
    T::AbstractArray{Float64}

    r::NTuple{N, AbstractArray{Float64}}
    k::NTuple{N, AbstractArray{Float64}}
end

function gen_grid_r(grid::Grid{N}) where N
    bounds, npoint = grid.bounds, grid.npoint

    axes = ntuple(N) do i
        Vector(range(bounds[i][1], step=(bounds[i][2] - bounds[i][1]) / npoint, length=npoint))
    end

    return ntuple(N) do i
        reshape(axes[i], ntuple(j -> j == i ? (:) : 1, N)...)
    end
end

function gen_grid_k(grid::Grid{N}) where N
    bounds, npoint = grid.bounds, grid.npoint

    axes = ntuple(N) do i
        Vector(fftfreq(npoint, 2pi * npoint / (bounds[i][2] - bounds[i][1])))
    end

    return ntuple(N) do i
        reshape(axes[i], ntuple(j -> j == i ? (:) : 1, N)...)
    end
end

function gen_wfn(ic::InitialConditions{N}, r::NTuple{N, AbstractArray{Float64}}, npoint::Int, nstate::Int) where N
    W = zeros(ComplexF64, (ntuple(_ -> npoint, N)..., nstate))

    components = ntuple(N) do i
        @. exp(-0.5 * ic.gamma[i] * (r[i] - ic.pos[i])^2 + im * ic.mom[i] * (r[i] - ic.pos[i]))
    end

    broadcast!(*, selectdim(W, N + 1, ic.state + 1), components...)

    return W
end

function get_dr(grid::Grid{N}) where N
    return prod(map(limit -> (limit[2] - limit[1]) / grid.npoint, grid.bounds))
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

function flat_prop(P::AbstractArray{ComplexF64})
    return reshape(P, size(P, 1), size(P, 2), :)
end

function calc_pe(W::AbstractArray{ComplexF64}, V::AbstractArray{Float64}, grid::Grid{N}) where N
    W_f, V_f = flat_wfn(W), flat_pot(V)

    @tullio pe := conj(W_f[i, I]) * V_f[I, J, i] * W_f[i, J]

    return real(pe) * get_dr(grid)
end

function calc_pops(W::AbstractArray{ComplexF64}, grid::Grid{N}) where N
    W_f = flat_wfn(W)
    
    @tullio pop[I] := abs2(W_f[i, I])
    
    return Tuple(pop .* get_dr(grid))
end

function normalize!(W::AbstractArray{ComplexF64}, grid::Grid{N}) where N
    W ./= calc_norm(W, grid) 
end

function calc_pos(W::AbstractArray{ComplexF64}, grid::Grid{N}, r::NTuple{N, AbstractArray{Float64}}) where N
    return ntuple(i -> sum(abs2.(W) .* r[i]) * get_dr(grid), N)
end

function to_kspace(W::AbstractArray{ComplexF64})
    return fft(W, 1:ndims(W) - 1)
end

function get_dk(grid::Grid{N}) where N
    return get_dr(grid) / grid.npoint^N
end

function calc_mom(W::AbstractArray{ComplexF64}, grid::Grid{N}, k::NTuple{N, AbstractArray{Float64}}) where N
    return ntuple(i -> sum(abs2.(W) .* k[i]) * get_dk(grid), N)
end

function calc_ke(W::AbstractArray{ComplexF64}, T::AbstractArray{Float64}, grid::Grid{N}) where N
    return sum(abs2.(W) .* T) * get_dk(grid)
end

function get_pot_eigen(V::AbstractArray{Float64, P}) where P
    A, U = similar(V, (size(V, 1), size(V)[3:end]...)), similar(V)

    for l in CartesianIndices(size(V)[3:end])
        A[:, l], U[:, :, l] = eigen(Symmetric(view(V, :, :, l)))
    end

    U_f = reshape(U, size(V, 1), size(V, 2), :)
    
    for i in 2:size(U_f, 3), j in 1:size(V, 1)
        if dot(view(U_f, :, j, i - 1), view(U_f, :, j, i)) < 0 U_f[:, j, i] .*= -1.0 end
    end
    
    return A, U
end

function get_prop_r(A::AbstractArray{Float64}, U::AbstractArray{Float64}, dt::ComplexF64)
    R = similar(U, ComplexF64)

    for l in CartesianIndices(size(A)[2:end])
        @views R[:, :, l] = U[:, :, l] * Diagonal(exp.(-im * A[:, l] * dt)) * U[:, :, l]'
    end

    return R
end

function get_prop_k(m::Float64, k::NTuple{N, AbstractArray{Float64}}, dt::ComplexF64) where N
    return exp.(-0.5im * reduce(.+, k[i].^2 / m for i in 1:N) * dt)
end

function apply_cap!(R::AbstractArray{ComplexF64}, absorber::Absorber{N}, r::NTuple{N, AbstractArray{Float64}}, dt::Float64) where N
    C = map(1:N) do i
        @. exp(absorber.exponent * max(0, absorber.bounds[i][1] - r[i], r[i] - absorber.bounds[i][2])) - 1
    end

    decay = exp.(-0.5 .* reduce(.+, C) .* dt)

    R .*= ComplexF64.(reshape(decay, 1, 1, size(decay)...))
end

function propagate_r!(W::AbstractArray{ComplexF64}, R::AbstractArray{ComplexF64})
    W_f, R_f = flat_wfn(W), flat_prop(R)

    @tullio W_t[i, I] := R_f[I, J, i] * W_f[i, J]

    W_f .= W_t
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

function propagate!(W::AbstractArray{ComplexF64}, R::AbstractArray{ComplexF64}, K::AbstractArray{ComplexF64}, U::AbstractArray{Float64}, grid::Grid{N}, adia::Bool) where N
    propagate_r!(W, R); d1 = calc_decay(W, R, U, grid, adia)

    propagate_k!(W, K)

    propagate_r!(W, R); d2 = calc_decay(W, R, U, grid, adia)

    return d1 .+ d2
end

function print_iter(i::Int, obs::Observables{N}, elapsed::UInt64) where N
    fmt_pos = join([@sprintf("%10.4f", x) for x in obs.pos], " ")
    fmt_mom = join([@sprintf("%10.4f", x) for x in obs.mom], " ")
    fmt_pop = join([@sprintf("%10.4f", x) for x in obs.pop], " ")

    values = (i - 1, obs.ke, obs.pe, obs.ke + obs.pe, fmt_pos, fmt_mom, fmt_pop, obs.norm, format_duration(elapsed))

    @printf("%7d %12.6f %12.6f %12.6f [%s] [%s] [%s] %9.4f %s\n", values...)
end

function print_header(state::Int, ndim::Int, nstate::Int, itp::Bool)
    dim_w, state_w = 11 * ndim + 1, 11 * nstate + 1

    itp && print("\nSTATE $state ITP")

    labels = ("ITER", "KIN (Eh)", "POT (Eh)", "TOT (Eh)", dim_w, "POS (a0)", dim_w, "MOM (hb/a0)", state_w, "POPULATION", "NORM", "TIME")

    @printf("\n%7s %12s %12s %12s %*s %*s %*s %9s %-s\n", labels...)
end

function format_duration(nanos::UInt64)
    return Dates.format(Time(0) + Nanosecond(nanos), "HH:MM:SS.sss")
end

function to_adia(W::AbstractArray{ComplexF64}, U::AbstractArray{Float64})
    W_f, U_f = flat_wfn(W), flat_pot(U)

    @tullio W_a[i, I] := conj(U_f[J, I, i]) * W_f[i, J]

    return reshape(W_a, size(W))
end

function to_dia(W::AbstractArray{ComplexF64}, U::AbstractArray{Float64})
    W_f, U_f = flat_wfn(W), flat_pot(U)

    @tullio W_d[i, I] := U_f[I, J, i] * W_f[i, J]

    return reshape(W_d, size(W))
end

function calc_observables(ctx::SimulationContext{N}, sim::Simulation, grid::Grid{N}, decay::Tuple, writer::Writer, log::Bool) where N
    need_k, need_a = log || !isnothing(writer.mom) || !isnothing(writer.ke) || !isnothing(writer.te), (log || !isnothing(writer.pop)) && sim.adia

    W_d, W_k, W_a = ctx.W, need_k ? to_kspace(ctx.W) : ctx.W, need_a ? to_adia(ctx.W, ctx.U) : ctx.W

    norm = log || !isnothing(writer.norm) ? calc_norm(W_d, grid)          : nothing
    pops = log || !isnothing(writer.pop)  ? calc_pops(W_a, grid) .+ decay : nothing

    pe = log || !isnothing(writer.pe) || !isnothing(writer.te) ? calc_pe(W_d, ctx.V, grid) : nothing
    ke = log || !isnothing(writer.ke) || !isnothing(writer.te) ? calc_ke(W_k, ctx.T, grid) : nothing

    pos = log || !isnothing(writer.pos) ? calc_pos(W_d, grid, ctx.r) : nothing
    mom = log || !isnothing(writer.mom) ? calc_mom(W_k, grid, ctx.k) : nothing

    return Observables{N, size(ctx.V, 1)}(pos=pos, mom=mom, pop=pops, norm=norm, pe=pe, ke=ke)
end

function init_history(writer::Writer, ndim::Int, nstate::Int)
    return History{ndim, nstate}(
        pos = !isnothing(writer.pos) ? NTuple{ndim,   Float64}[] : nothing,
        mom = !isnothing(writer.mom) ? NTuple{ndim,   Float64}[] : nothing,
        pop = !isnothing(writer.pop) ? NTuple{nstate, Float64}[] : nothing,
        
        norm = !isnothing(writer.norm) ? Float64[] : nothing,
        pe   = !isnothing(writer.pe)   ? Float64[] : nothing,
        ke   = !isnothing(writer.ke)   ? Float64[] : nothing,

        wfn = !isnothing(writer.wfn) ? AbstractArray{ComplexF64}[] : nothing
    )
end

function get_kin_op(m::Float64, k::NTuple{N, AbstractArray{Float64}}) where N
    return 0.5 * reduce(.+, k[i].^2 for i in 1:N) / m
end

function overlap(W1::AbstractArray{ComplexF64}, W2::AbstractArray{ComplexF64}, grid::Grid{N}) where N
    return sum(conj(W1[l]) * W2[l] for l in CartesianIndices(W1)) * get_dr(grid)
end

function project_out!(W1::AbstractArray{ComplexF64}, W2::AbstractArray{ComplexF64}, grid::Grid{N}) where N
    W1 .-= overlap(W2, W1, grid) .* W2
end

function log_final_pop(W::AbstractArray{ComplexF64}, U::AbstractArray{Float64}, grid::Grid{N}, decay::Tuple, adia::Bool) where N
    pop = calc_pops(adia ? to_adia(W, U) : W, grid) .+ decay

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

function parse_config(config::Dict{String, Any})
    dt::Float64 = config["simulation"]["time_step" ]
    iters::Int  = config["simulation"]["iterations"]

    adia::Bool        = get(config["simulation"], "adiabatic",    false)
    itp::Int          = get(config["simulation"], "imaginary",    0    )
    log_interval::Int = get(config["simulation"], "log_interval", 1    )
    m::Float64        = get(config["simulation"], "mass",         1    )

    sim = Simulation(itp=itp, iters=iters, dt=dt, m=m, adia=adia, log_interval=log_interval)

    bounds = Tuple(Tuple{Float64, Float64}(b) for b in config["grid"]["bounds"])

    grid = Grid(bounds=bounds, npoint=config["grid"]["npoint"])

    pos   = Tuple(Float64(q) for q in config["initial_conditions"]["position"])
    mom   = Tuple(Float64(p) for p in config["initial_conditions"]["momentum"])
    gamma = Tuple(Float64(g) for g in config["initial_conditions"]["gamma"   ])

    state::Int    = get(config["initial_conditions"], "state",     0    )
    ic_adia::Bool = get(config["initial_conditions"], "adiabatic", false)

    ic = InitialConditions(pos=pos, mom=mom, gamma=gamma, state=state, adia=ic_adia)

    writer = if haskey(config, "write")
        pos  = get(config["write"], "position",         nothing)
        mom  = get(config["write"], "momentum",         nothing)
        pop  = get(config["write"], "population",       nothing)
        pe   = get(config["write"], "potential_energy", nothing)
        ke   = get(config["write"], "kinetic_energy",   nothing)
        te   = get(config["write"], "total_energy",     nothing)
        norm = get(config["write"], "norm",             nothing)
        wfn  = get(config["write"], "wavefunction",     nothing)

        Writer(pos=pos, mom=mom, pop=pop, pe=pe, ke=ke, te=te, norm=norm, wfn=wfn)
    else Writer() end

    absorber = if haskey(config, "cap")
        bounds = Tuple(Tuple{Float64, Float64}(b) for b in config["cap"]["bounds"])

        Absorber(bounds=bounds, exponent=config["cap"]["exponent"])
    else nothing end

    return sim, grid, ic, POTENTIALS[config["potential"]["name"]], absorber, writer
end

function log_history!(history::History{N, S}, W::AbstractArray{ComplexF64}, obs::Observables{N, S}) where {N, S}
    !isnothing(history.pos) && push!(history.pos, obs.pos)
    !isnothing(history.mom) && push!(history.mom, obs.mom)
    !isnothing(history.pop) && push!(history.pop, obs.pop)
    
    !isnothing(history.norm) && push!(history.norm, obs.norm)
    !isnothing(history.pe  ) && push!(history.pe,   obs.pe  )
    !isnothing(history.ke  ) && push!(history.ke,   obs.ke  )

    !isnothing(history.wfn) && push!(history.wfn, copy(W))
end

function export_history(history::History{N, S}, ctx::SimulationContext{N}, sim::Simulation, writer::Writer, state::Int) where {N, S}
    to_matrix(v::Vector{NTuple{M, Float64}}) where M = [row[col] for row in v, col in 1:M]

    function format_name(fname::String, state::Int, itp::Bool)
        return itp ? "$(splitext(fname)[1])_STATE-$(state)$(splitext(fname)[2])" : fname
    end

    t = collect(0:sim.dt:sim.dt * sim.iters)

    !isnothing(history.pos) && write_matrix(format_name(writer.pos, state, sim.itp > 0), hcat(t, to_matrix(history.pos)))
    !isnothing(history.mom) && write_matrix(format_name(writer.mom, state, sim.itp > 0), hcat(t, to_matrix(history.mom)))
    !isnothing(history.pop) && write_matrix(format_name(writer.pop, state, sim.itp > 0), hcat(t, to_matrix(history.pop)))

    !isnothing(history.norm) && write_matrix(format_name(writer.norm, state, sim.itp > 0), hcat(t, history.norm))

    !isnothing(history.pe) && write_matrix(format_name(writer.pe, state, sim.itp > 0), hcat(t, history.pe              ))
    !isnothing(history.ke) && write_matrix(format_name(writer.ke, state, sim.itp > 0), hcat(t, history.ke              ))
    !isnothing(writer.te ) && write_matrix(format_name(writer.te, state, sim.itp > 0), hcat(t, history.pe .+ history.ke))

    if !isnothing(history.wfn)
        write_matrix(format_name(writer.wfn, state, sim.itp > 0), hcat(flatten_coords(ctx.r), flatten_wfn_history(history.wfn)))
    end
end

function flatten_coords(r::NTuple{N, AbstractArray{Float64}}) where N
    return [map(vec, r)[d][c[d]] for c in vec(CartesianIndices(map(length, map(vec, r)))), d in 1:N]
end

function flatten_wfn_history(wfns::Vector{<:AbstractArray{ComplexF64}})
    reduce(hcat, [hcat(real.(vec(slice)), imag.(vec(slice))) for W in wfns for slice in eachslice(W, dims=ndims(W))])
end

function write_matrix(fname::String, mat::AbstractArray{Float64})
    h5open(fname, "w") do file file["data"] = mat end
end

function calc_decay(W::AbstractArray{ComplexF64}, R::AbstractArray{ComplexF64}, U::AbstractArray{Float64}, grid::Grid{N}, adia::Bool) where N
    W_f, R_f = flat_wfn(adia ? to_adia(W, U) : W), flat_prop(R)

    @tullio surv[i] := abs2(R_f[I, 1, i])

    @tullio decay[I] := abs2(W_f[i, I]) * (1 / max(surv[i], 1e-14) - 1)

    return Tuple(decay .* get_dr(grid))
end

function run_qd(config::Dict{String, Any}, enable_print::Bool = true)
    sim, grid, ic, pot, absorber, writer = parse_config(config)

    @timeit "INITIALIZATION" begin
        r = gen_grid_r(grid)
        k = gen_grid_k(grid)

        V, m = pot.fn(r...), sim.m; T = get_kin_op(m, k)

        A, U = get_pot_eigen(V)

        R = get_prop_r(A, U, 0.5 * (sim.itp > 0 ? -im * sim.dt : sim.dt + 0im))
        K = get_prop_k(m, k, 1.0 * (sim.itp > 0 ? -im * sim.dt : sim.dt + 0im))

        if !isnothing(absorber)
            apply_cap!(R, absorber, r, 0.5 * sim.dt)
        end
    end

    opt_wfn, opt_wfn_te, output = AbstractArray{ComplexF64}[], Float64[], Observables{length(r), pot.nstate}[]

    for i in 1:(sim.itp > 0 ? sim.itp : 1)
        W = gen_wfn(ic, r, grid.npoint, pot.nstate); normalize!(W, grid)

        if ic.adia
            W = to_dia(W, U)
        end

        ctx, decay = SimulationContext(W=W, R=R, K=K, V=V, A=A, U=U, T=T, r=r, k=k), ntuple(_ -> 0.0, pot.nstate)

        enable_print && print_header(i, length(r), pot.nstate, sim.itp > 0)

        history, start_time = init_history(writer, length(r), pot.nstate), time_ns()

        for j in 1:sim.iters + 1
            log = (j - 1) % sim.log_interval == 0 || j == sim.iters + 1

            if j > 1
                @timeit "PROPAGATION" if isnothing(absorber) propagate!(W, R, K) else decay = decay .+ propagate!(W, R, K, U, grid, sim.adia) end
            end

            if sim.itp > 0
                @timeit "PROJECTION" project_and_normalize!(W, opt_wfn, grid)
            end

            if log || writer != Writer()
                @timeit "OBSERVABLES" obs = calc_observables(ctx, sim, grid, decay, writer, log)

                enable_print && log && print_iter(j, obs, time_ns() - start_time)

                if j == sim.iters + 1
                    sim.itp > 1 && push!(opt_wfn_te, obs.pe + obs.ke); push!(output, obs)
                end

                if log start_time = time_ns() end

                if writer != Writer()
                    @timeit "HISTORY APPEND" log_history!(history, sim.adia ? to_adia(W, U) : W, obs)
                end
            end
        end

        if pot.nstate > 1 && sim.itp == 0
            enable_print && log_final_pop(W, ctx.U, grid, decay, sim.adia)
        end

        sim.itp > 0 && push!(opt_wfn, W)

        if writer != Writer()
            @timeit "DATA EXPORT" export_history(history, ctx, sim, writer, i)
        end
    end

    sim.itp > 1 && enable_print && log_final_te(opt_wfn_te)

    return output
end

end # module QuantumDynamics
