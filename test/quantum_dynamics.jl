using Test, Zinq.QuantumDynamics

@testset "Tully-1 (Real Time)" begin
    config = Dict{String, Any}(
        "simulation" => Dict{String, Any}(
            "iterations" => 3000,
            "time_step" => 1,
            "mass" => 2000.0,
            "adiabatic" => true,
            "log_interval" => 500
        ),
        "grid" => Dict{String, Any}(
            "bounds" => [[-24.0, 32.0]],
            "npoint" => 512
        ),
        "initial_conditions" => Dict{String, Any}(
            "position" => [-10],
            "momentum" => [15],
            "gamma" => [2],
            "state" => 1,
            "adiabatic" => true
        ),
        "potential" => Dict{String, Any}(
            "name" => "tully_1"
        )
    )

    output = run_qd(config, false)

    @test isapprox(output[1].pos[1], 13.41124587399755, atol=1e-12)
    @test isapprox(output[1].mom[1], 16.01019634771595, atol=1e-12)
    @test isapprox(output[1].pop[1],  0.41038519680951, atol=1e-12)
    @test isapprox(output[1].pop[2],  0.58961480319084, atol=1e-12)
    @test isapprox(output[1].norm,    1.00000000000017, atol=1e-12)
    @test isapprox(output[1].pe,      0.00179229580696, atol=1e-12)
    @test isapprox(output[1].ke,      0.06470770264329, atol=1e-12)
end
