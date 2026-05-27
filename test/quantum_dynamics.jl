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
            "npoint" => 2048
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

    @test isapprox(output[1].pos[1], 13.41125203921865, atol=1e-12)
    @test isapprox(output[1].mom[1], 16.01020084897931, atol=1e-12)
    @test isapprox(output[1].pop[1],  0.41038326889215, atol=1e-12)
    @test isapprox(output[1].pop[2],  0.58961673110844, atol=1e-12)
    @test isapprox(output[1].norm,    1.00000000000030, atol=1e-12)
    @test isapprox(output[1].pe,      0.00179233436531, atol=1e-12)
    @test isapprox(output[1].ke,      0.06470766408496, atol=1e-12)
end

@testset "Harmonic-1D (Imaginary Time)" begin
    config = Dict{String, Any}(
        "simulation" => Dict{String, Any}(
            "imaginary" => 2,
            "iterations" => 1000,
            "time_step" => 0.01,
            "mass" => 1,
        ),
        "grid" => Dict{String, Any}(
            "bounds" => [[-8.0, 8.0]],
            "npoint" => 256,
        ),
        "initial_conditions" => Dict{String, Any}(
            "position" => [1],
            "momentum" => [0],
            "gamma" => [2],
        ),
        "potential" => Dict{String, Any}(
            "name" => "harmonic",
        )
    )

    output = run_qd(config, false)

    @test isapprox(output[1].pos[1], 0.00006053550965, atol=1e-12)
    @test isapprox(output[1].mom[1], 0.00000000000000, atol=1e-12)
    @test isapprox(output[1].norm,   1.00000000000000, atol=1e-12)
    @test isapprox(output[1].pe,     0.24999687654732, atol=1e-12)
    @test isapprox(output[1].ke,     0.25000312532402, atol=1e-12)

    @test isapprox(output[2].pos[1], -0.00002269862596, atol=1e-12)
    @test isapprox(output[2].mom[1],  0.00000000000000, atol=1e-12)
    @test isapprox(output[2].norm,    1.00000000000000, atol=1e-12)
    @test isapprox(output[2].pe,      0.74999062375148, atol=1e-12)
    @test isapprox(output[2].ke,      0.75000937489134, atol=1e-12)
end

@testset "Harmonic-1D (Real Time)" begin
    config = Dict{String, Any}(
        "simulation" => Dict{String, Any}(
            "iterations" => 1000,
            "time_step" => 0.01,
            "mass" => 1,
        ),
        "grid" => Dict{String, Any}(
            "bounds" => [[-8.0, 8.0]],
            "npoint" => 256,
        ),
        "initial_conditions" => Dict{String, Any}(
            "position" => [1],
            "momentum" => [0],
            "gamma" => [2],
        ),
        "potential" => Dict{String, Any}(
            "name" => "harmonic",
        )
    )

    output = run_qd(config, false)

    @test isapprox(output[1].pos[1], -0.83904886054688, atol=1e-12)
    @test isapprox(output[1].mom[1],  0.54404927138077, atol=1e-12)
    @test isapprox(output[1].norm,    1.00000000000007, atol=1e-12)
    @test isapprox(output[1].pe,      0.58800407385449, atol=1e-12)
    @test isapprox(output[1].ke,      0.53699500124751, atol=1e-12)
end

@testset "Harmonic-2D (Imaginary Time)" begin
    config = Dict{String, Any}(
        "simulation" => Dict{String, Any}(
            "imaginary" => 2,
            "iterations" => 1000,
            "time_step" => 0.01,
            "mass" => 1,
        ),
        "grid" => Dict{String, Any}(
            "bounds" => [[-8.0, 8.0], [-8.0, 8.0]],
            "npoint" => 128,
        ),
        "initial_conditions" => Dict{String, Any}(
            "position" => [1, 1],
            "momentum" => [0, 0],
            "gamma" => [2, 2],
        ),
        "potential" => Dict{String, Any}(
            "name" => "harmonic",
        )
    )

    output = run_qd(config, false)

    @test isapprox(output[1].pos[1], 0.00006053550965, atol=1e-12)
    @test isapprox(output[1].pos[2], 0.00006053550965, atol=1e-12)
    @test isapprox(output[1].mom[1], 0.00000000000000, atol=1e-12)
    @test isapprox(output[1].mom[2], 0.00000000000000, atol=1e-12)
    @test isapprox(output[1].norm,   1.00000000000000, atol=1e-12)
    @test isapprox(output[1].pe,     0.49999375309465, atol=1e-12)
    @test isapprox(output[1].ke,     0.50000625064804, atol=1e-12)

    @test isapprox(output[2].pos[1], -0.00001134793993, atol=1e-12)
    @test isapprox(output[2].pos[2], -0.00001134793993, atol=1e-12)
    @test isapprox(output[2].mom[1],  0.00000000000000, atol=1e-12)
    @test isapprox(output[2].mom[2],  0.00000000000000, atol=1e-12)
    @test isapprox(output[2].norm,    1.00000000000000, atol=1e-12)
    @test isapprox(output[2].pe,      0.99998749835204, atol=1e-12)
    @test isapprox(output[2].ke,      1.00001249941377, atol=1e-12)
end

@testset "Harmonic-2D (Real Time)" begin
    config = Dict{String, Any}(
        "simulation" => Dict{String, Any}(
            "iterations" => 1000,
            "time_step" => 0.01,
            "mass" => 1,
        ),
        "grid" => Dict{String, Any}(
            "bounds" => [[-8.0, 8.0], [-8.0, 8.0]],
            "npoint" => 128,
        ),
        "initial_conditions" => Dict{String, Any}(
            "position" => [1, 1],
            "momentum" => [0, 0],
            "gamma" => [2, 2],
        ),
        "potential" => Dict{String, Any}(
            "name" => "harmonic",
        )
    )

    output = run_qd(config, false)

    @test isapprox(output[1].pos[1], -0.83904886054684, atol=1e-12)
    @test isapprox(output[1].pos[2], -0.83904886054684, atol=1e-12)
    @test isapprox(output[1].mom[1],  0.54404927138074, atol=1e-12)
    @test isapprox(output[1].mom[2],  0.54404927138074, atol=1e-12)
    @test isapprox(output[1].norm,    1.00000000000004, atol=1e-12)
    @test isapprox(output[1].pe,      1.17600814770893, atol=1e-12)
    @test isapprox(output[1].ke,      1.07399000249496, atol=1e-12)
end
