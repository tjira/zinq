using Test, Zinq

@testset "Zinq.jl" begin
    @testset "Quantum Dynamics" begin
        include("quantum_dynamics.jl")
    end
end
