module QuantumDynamics

include("Potentials.jl")

using .Potentials, FFTW, HDF5, TOML, Dates, LinearAlgebra, Printf, TimerOutputs, Tullio

export parse_config_qd, run_qd, run_qd_stable

Base.@kwdef struct Absorber{N}
    bounds::NTuple{N, Tuple{Float64, Float64}}; exponent::Float64
end

Base.@kwdef struct Grid{N, T <: AbstractArray{Float64}}
    bounds::NTuple{N, Tuple{Float64, Float64}}; npoint::Int

    r::NTuple{N, T}
    k::NTuple{N, T}
end

Base.@kwdef mutable struct Wavefunction{T <: AbstractArray{ComplexF64}}
    data::T; decay::Vector{Float64}
end

Base.@kwdef struct Hamiltonian{TV <: AbstractArray{Float64}, TK <: AbstractArray{Float64}, TA <: AbstractArray{Float64}}
    V::TV; K::TK; A::TA; U::TV
end

Base.@kwdef struct Propagator{TR <: AbstractArray{ComplexF64}, TK <: AbstractArray{ComplexF64}}
    R::TR; K::TK
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

Base.@kwdef struct History{N, S, T <: AbstractArray{ComplexF64}}
    pos::Union{Vector{NTuple{N, Float64}}, Nothing} = nothing
    mom::Union{Vector{NTuple{N, Float64}}, Nothing} = nothing
    pop::Union{Vector{NTuple{S, Float64}}, Nothing} = nothing

    norm::Union{Vector{Float64}, Nothing} = nothing
    pe  ::Union{Vector{Float64}, Nothing} = nothing
    ke  ::Union{Vector{Float64}, Nothing} = nothing

    wfn::Union{Vector{T}, Nothing} = nothing
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

function gen_grid_r(bounds::NTuple{N, Tuple{Float64, Float64}}, npoint::Int) where N
    axes = ntuple(N) do i
        Vector(range(bounds[i][1], step=(bounds[i][2] - bounds[i][1]) / npoint, length=npoint))
    end

    return ntuple(N) do i
        reshape(axes[i], ntuple(j -> j == i ? (:) : 1, N)...)
    end
end

function gen_grid_k(bounds::NTuple{N, Tuple{Float64, Float64}}, npoint::Int) where N
    axes = ntuple(N) do i
        Vector(fftfreq(npoint, 2pi * npoint / (bounds[i][2] - bounds[i][1])))
    end

    return ntuple(N) do i
        reshape(axes[i], ntuple(j -> j == i ? (:) : 1, N)...)
    end
end

function Grid(bounds::NTuple{N, Tuple{Float64, Float64}}, npoint::Int) where N
    r = gen_grid_r(bounds, npoint)
    k = gen_grid_k(bounds, npoint)

    return Grid(bounds=bounds, npoint=npoint, r=r, k=k)
end

function gen_wfn(ic::InitialConditions{N}, grid::Grid{N}, nstate::Int) where N
    W = zeros(ComplexF64, (ntuple(_ -> grid.npoint, N)..., nstate))

    components = ntuple(N) do i
        @. exp(-0.5 * ic.gamma[i] * (grid.r[i] - ic.pos[i])^2 + im * ic.mom[i] * (grid.r[i] - ic.pos[i]))
    end

    broadcast!(*, selectdim(W, N + 1, ic.state + 1), components...)

    return Wavefunction(data=W, decay=zeros(Float64, nstate))
end

function get_dr(grid::Grid{N}) where N
    return prod(map(limit -> (limit[2] - limit[1]) / grid.npoint, grid.bounds))
end

function calc_norm(psi::Wavefunction, grid::Grid{N}) where N
    return sqrt(sum(abs2, psi.data) * get_dr(grid))
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

function calc_pe(psi::Wavefunction, ham::Hamiltonian, grid::Grid{N}) where N
    W_f, V_f = flat_wfn(psi.data), flat_pot(ham.V)

    @tullio pe := conj(W_f[i, I]) * V_f[I, J, i] * W_f[i, J]

    return real(pe) * get_dr(grid)
end

function calc_pops(psi::Wavefunction, grid::Grid{N}) where N
    W_f = flat_wfn(psi.data)
    
    @tullio pop[I] := abs2(W_f[i, I])
    
    return Vector(pop .* get_dr(grid))
end

function normalize!(psi::Wavefunction{T}, grid::Grid{N}) where {N, T}
    psi.data ./= calc_norm(psi, grid) 
end

function calc_pos(psi::Wavefunction, grid::Grid{N}) where N
    return ntuple(i -> sum(abs2.(psi.data) .* grid.r[i]) * get_dr(grid), N)
end

function to_kspace(psi::Wavefunction)
    return Wavefunction(data=fft(psi.data, 1:ndims(psi.data) - 1), decay=psi.decay)
end

function get_dk(grid::Grid{N}) where N
    return get_dr(grid) / grid.npoint^N
end

function calc_mom(psi_k::Wavefunction, grid::Grid{N}) where N
    return ntuple(i -> sum(abs2.(psi_k.data) .* grid.k[i]) * get_dk(grid), N)
end

function calc_ke(psi_k::Wavefunction, ham::Hamiltonian, grid::Grid{N}) where N
    return sum(abs2.(psi_k.data) .* ham.K) * get_dk(grid)
end

function get_pot_eigen(V::T) where T
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

function get_prop_k(m::Float64, grid::Grid{N}, dt::ComplexF64) where N
    return exp.(-0.5im * reduce(.+, grid.k[i].^2 / m for i in 1:N) * dt)
end

function apply_cap!(prop::Propagator, absorber::Absorber{N}, grid::Grid{N}, dt::Float64) where N
    C = map(1:N) do i
        @. exp(absorber.exponent * max(0, absorber.bounds[i][1] - grid.r[i], grid.r[i] - absorber.bounds[i][2])) - 1
    end

    decay = exp.(-0.5 .* reduce(.+, C) .* dt)

    prop.R .*= ComplexF64.(reshape(decay, 1, 1, size(decay)...))
end

function propagate_r!(psi::Wavefunction, prop::Propagator)
    W_f, R_f = flat_wfn(psi.data), flat_prop(prop.R)

    @tullio W_t[i, I] := R_f[I, J, i] * W_f[i, J]

    W_f .= W_t
end

function propagate_k!(psi::Wavefunction, prop::Propagator)
    fft!(psi.data, 1:ndims(psi.data) - 1)

    psi.data .*= prop.K

    ifft!(psi.data, 1:ndims(psi.data) - 1)
end

function propagate!(psi::Wavefunction, prop::Propagator, ham::Hamiltonian, grid::Grid{N}, adia::Bool, absorber = nothing) where N
    propagate_r!(psi, prop)

    if !isnothing(absorber)
        psi.decay .+= calc_decay(psi, prop, ham, grid, adia)
    end

    propagate_k!(psi, prop)
    propagate_r!(psi, prop)

    if !isnothing(absorber)
        psi.decay .+= calc_decay(psi, prop, ham, grid, adia)
    end
end

function print_iter(i::Int, obs::Observables{N, S}, elapsed::UInt64) where {N, S}
    fmt_pos = join([@sprintf("%10.4f", x) for x in obs.pos], " ")
    fmt_mom = join([@sprintf("%10.4f", x) for x in obs.mom], " ")
    fmt_pop = join([@sprintf("%10.4f", x) for x in obs.pop], " ")

    values = (i - 1, obs.ke, obs.pe, obs.ke + obs.pe, fmt_pos, fmt_mom, fmt_pop, obs.norm, format_duration(elapsed))

    @printf("%7d %12.6f %12.6f %12.6f [%s] [%s] [%s] %9.4f %s\n", values...)
end

function print_header(state::Int, itp::Bool, ::Val{N}, ::Val{S}) where {N, S}
    dim_w, state_w = 11 * N + 1, 11 * S + 1

    itp && print("\nSTATE $state ITP")

    labels = ("ITER", "KIN (Eh)", "POT (Eh)", "TOT (Eh)", dim_w, "POS (a0)", dim_w, "MOM (hb/a0)", state_w, "POPULATION", "NORM", "TIME")

    @printf("\n%7s %12s %12s %12s %*s %*s %*s %9s %-s\n", labels...)
end

function format_duration(nanos::UInt64)
    return Dates.format(Time(0) + Nanosecond(nanos), "HH:MM:SS.sss")
end

function to_adia(psi::Wavefunction, ham::Hamiltonian)
    W_f, U_f = flat_wfn(psi.data), flat_pot(ham.U)

    @tullio W_a[i, I] := conj(U_f[J, I, i]) * W_f[i, J]

    return Wavefunction(data=reshape(W_a, size(psi.data)), decay=psi.decay)
end

function to_dia(psi::Wavefunction, ham::Hamiltonian)
    W_f, U_f = flat_wfn(psi.data), flat_pot(ham.U)

    @tullio W_d[i, I] := U_f[I, J, i] * W_f[i, J]

    return Wavefunction(data=reshape(W_d, size(psi.data)), decay=psi.decay)
end

function calc_observables(psi::Wavefunction, ham::Hamiltonian, grid::Grid{N}, sim::Simulation, writer::Writer, log::Bool) where N
    need_k, need_a = log || !isnothing(writer.mom) || !isnothing(writer.ke) || !isnothing(writer.te), (log || !isnothing(writer.pop)) && sim.adia

    psi_k, psi_a = need_k ? to_kspace(psi) : nothing, need_a ? to_adia(psi, ham) : nothing

    pops = log || !isnothing(writer.pop) ? Tuple(calc_pops(need_a ? psi_a : psi, grid) .+ psi.decay) : nothing

    pe = log || !isnothing(writer.pe) || !isnothing(writer.te) ? calc_pe(psi,   ham, grid) : nothing
    ke = log || !isnothing(writer.ke) || !isnothing(writer.te) ? calc_ke(psi_k, ham, grid) : nothing

    norm = log || !isnothing(writer.norm) ? calc_norm(psi,   grid) : nothing
    pos  = log || !isnothing(writer.pos ) ? calc_pos( psi,   grid) : nothing
    mom  = log || !isnothing(writer.mom ) ? calc_mom( psi_k, grid) : nothing

    return Observables{N, size(ham.V, 1)}(pos=pos, mom=mom, pop=pops, norm=norm, pe=pe, ke=ke)
end

function init_history(writer::Writer, ::Val{N}, ::Val{S}, ::Type{T}) where {N, S, T}
    return History{N, S, T}(
        pos = !isnothing(writer.pos) ? NTuple{N, Float64}[] : nothing,
        mom = !isnothing(writer.mom) ? NTuple{N, Float64}[] : nothing,
        pop = !isnothing(writer.pop) ? NTuple{S, Float64}[] : nothing,
        
        norm = !isnothing(writer.norm) ? Float64[] : nothing,
        pe   = !isnothing(writer.pe)   ? Float64[] : nothing,
        ke   = !isnothing(writer.ke)   ? Float64[] : nothing,

        wfn = !isnothing(writer.wfn) ? T[] : nothing
    )
end

function get_kin_op(m::Float64, grid::Grid{N}) where N
    return 0.5 * reduce(.+, grid.k[i].^2 for i in 1:N) / m
end

function overlap(psi1::Wavefunction, psi2::Wavefunction, grid::Grid{N}) where N
    return sum(conj(psi1.data[l]) * psi2.data[l] for l in CartesianIndices(psi1.data)) * get_dr(grid)
end

function project_out!(psi1::Wavefunction, psi2::Wavefunction, grid::Grid{N}) where N
    psi1.data .-= overlap(psi2, psi1, grid) .* psi2.data
end

function log_final_pop(psi::Wavefunction, ham::Hamiltonian, grid::Grid{N}, adia::Bool) where N
    pop = calc_pops(adia ? to_adia(psi, ham) : psi, grid) .+ psi.decay

    for (i, p) in enumerate(pop)
        @printf("%sFINAL POPULATION OF STATE %02d: %.6f\n", i == 1 ? "\n" : "", i, p)
    end
end

function log_final_te(tes::Vector{Float64})
    for (i, te) in enumerate(tes)
        @printf("%sENERGY OF STATE %02d: %8.6f Eh\n", i == 1 ? "\n" : "", i, te)
    end
end

function project_and_normalize!(psi::Wavefunction, wfns::Vector{Wavefunction}, grid::Grid{N}) where N
    for psi2 in wfns
        project_out!(psi, psi2, grid)
    end

    normalize!(psi, grid)
end

function log_history!(history::History{N, S}, psi::Wavefunction, obs::Observables{N, S}) where {N, S}
    !isnothing(history.pos) && push!(history.pos, obs.pos)
    !isnothing(history.mom) && push!(history.mom, obs.mom)
    !isnothing(history.pop) && push!(history.pop, obs.pop)
    
    !isnothing(history.norm) && push!(history.norm, obs.norm)
    !isnothing(history.pe  ) && push!(history.pe,   obs.pe  )
    !isnothing(history.ke  ) && push!(history.ke,   obs.ke  )

    !isnothing(history.wfn) && push!(history.wfn, copy(psi.data))
end

function export_history(history::History{N, S}, grid::Grid{N}, sim::Simulation, writer::Writer, state::Int) where {N, S}
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
        write_matrix(format_name(writer.wfn, state, sim.itp > 0), hcat(flatten_coords(grid.r), flatten_wfn_history(history.wfn)))
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

function calc_decay(psi::Wavefunction, prop::Propagator, ham::Hamiltonian, grid::Grid{N}, adia::Bool) where N
    W_f, R_f = flat_wfn(adia ? to_adia(psi, ham).data : psi.data), flat_prop(prop.R)

    @tullio surv[i] := abs2(R_f[I, 1, i])

    @tullio decay[I] := abs2(W_f[i, I]) * (1 / max(surv[i], 1e-14) - 1)

    return Vector(decay .* get_dr(grid))
end

Base.@kwdef struct StableParams{N, S, T <: AbstractArray{Float64}, F <: Function}
    sim      ::Simulation
    grid     ::Grid{N, T}
    ic       ::InitialConditions{N}
    pot      ::Potential{S, F}
    absorber ::Union{Absorber{N}, Nothing}
    writer   ::Writer
end

function parse_config_qd(config::Dict{String, Any})
    dt::Float64 = config["simulation"]["time_step" ]
    iters::Int  = config["simulation"]["iterations"]

    adia::Bool        = get(config["simulation"], "adiabatic",    false)
    itp::Int          = get(config["simulation"], "imaginary",    0    )
    log_interval::Int = get(config["simulation"], "log_interval", 1    )
    m::Float64        = get(config["simulation"], "mass",         1    )

    sim = Simulation(itp=itp, iters=iters, dt=dt, m=m, adia=adia, log_interval=log_interval)

    bounds = Tuple(Tuple{Float64, Float64}(b) for b in config["grid"]["bounds"])

    grid = Grid(bounds, config["grid"]["npoint"])

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

    return StableParams(sim=sim, grid=grid, ic=ic, pot=POTENTIALS[config["potential"]["name"]], absorber=absorber, writer=writer)
end

function parse_config_qd(fname::String)
    return parse_config_qd(TOML.parsefile(fname))
end

function run_qd(config::Dict{String, Any}, enable_print::Bool = true)
    return run_qd_stable(parse_config_qd(config), enable_print)
end

function run_qd_stable(sp::StableParams{N, S, T, F}, enable_print::Bool = true) where {N, S, T, F}
    (;sim, grid, ic, pot, absorber, writer) = sp

    @timeit "INITIALIZATION" begin
        V, K = pot.fn(grid.r...), get_kin_op(sim.m, grid)

        A, U = get_pot_eigen(V); ham = Hamiltonian(V=V, K=K, A=A, U=U)

        R = get_prop_r(A,     U,    0.5 * (sim.itp > 0 ? -im * sim.dt : sim.dt + 0im))
        K = get_prop_k(sim.m, grid, 1.0 * (sim.itp > 0 ? -im * sim.dt : sim.dt + 0im))

        prop = Propagator(R=R, K=K)

        if !isnothing(absorber)
            apply_cap!(prop, absorber, grid, 0.5 * sim.dt)
        end
    end

    opt_psi, opt_psi_te, output = Wavefunction[], Float64[], Observables{N, S}[]

    for i in 1:(sim.itp > 0 ? sim.itp : 1)
        psi = gen_wfn(ic, grid, S); normalize!(psi, grid)

        if ic.adia
            psi = to_dia(psi, ham)
        end

        enable_print && print_header(i, sim.itp > 0, Val(N), Val(S))

        history, start_time = init_history(writer, Val(N), Val(S), typeof(psi.data)), time_ns()

        for j in 1:sim.iters + 1
            log = (j - 1) % sim.log_interval == 0 || j == sim.iters + 1

            if j > 1
                @timeit "PROPAGATION" propagate!(psi, prop, ham, grid, sim.adia, absorber)
            end

            if sim.itp > 0
                @timeit "PROJECTION" project_and_normalize!(psi, opt_psi, grid)
            end

            if log || writer != Writer()
                @timeit "OBSERVABLES" obs = calc_observables(psi, ham, grid, sim, writer, log)

                enable_print && log && print_iter(j, obs, time_ns() - start_time)

                if j == sim.iters + 1
                    sim.itp > 1 && push!(opt_psi_te, obs.pe + obs.ke); push!(output, obs)
                end

                if log start_time = time_ns() end

                if writer != Writer()
                    @timeit "HISTORY APPEND" log_history!(history, sim.adia ? to_adia(psi, ham) : psi, obs)
                end
            end
        end

        if S > 1 && sim.itp == 0
            enable_print && log_final_pop(psi, ham, grid, sim.adia)
        end

        sim.itp > 0 && push!(opt_psi, psi)

        if writer != Writer()
            @timeit "DATA EXPORT" export_history(history, grid, sim, writer, i)
        end
    end

    sim.itp > 1 && enable_print && log_final_te(opt_psi_te)

    return output
end

end # module QuantumDynamics
