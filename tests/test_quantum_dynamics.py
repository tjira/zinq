import pytest

from zinq.quantum_dynamics import Options, run


def error_message(observable: str, expected, actual) -> str:
    return f"INCORRECT {observable}: EXPECTED {expected}, GOT {actual}"


def test_quantum_dynamics_tully_1_real():
    options = {
        "type": "quantum_dynamics",
        "grid": {
            "limits": [[-24, 24]],
            "npoint": 512
        },
        "initial_conditions": {
            "momentum": [15],
            "position": [-10],
            "gamma": [2],
            "state": 1
        },
        "hamiltonian" : {
            "potential": {
                "name": "tully_1"
            },
            "mass": 2000,
        },
        "iterations": 3000,
        "time_step": 1,
        "log_interval": 500
    }

    result = run(Options.model_validate(options))

    expected_total_energies     = [0.06649999845026]
    expected_kinetic_energies   = [0.06470768763626]
    expected_potential_energies = [0.00179231081400]
    expected_norm_values        = [1.00000000000053]

    expected_position_values = [
        [13.41125215693357]
    ]

    expected_momentum_values = [
        [16.01020146127416]
    ]

    expected_populations = [
        [0.58961548497405, 0.41038451502648]
    ]

    for i, state in enumerate(result.states):
        actual = {
            "total_energy": state.total_energy,
            "kinetic_energy": state.kinetic_energy,
            "potential_energy": state.potential_energy,
            "position": state.position.tolist(),
            "momentum": state.momentum.tolist(),
            "norm": state.norm,
            "population": state.population.tolist(),
        }
        
        expected = {
            "total_energy": expected_total_energies[i],
            "kinetic_energy": expected_kinetic_energies[i],
            "potential_energy": expected_potential_energies[i],
            "position": expected_position_values[i],
            "momentum": expected_momentum_values[i],
            "norm": expected_norm_values[i],
            "population": expected_populations[i],
        }
        
        conditions = [actual[key] == pytest.approx(expected[key], abs=1e-08) for key in actual]

        for j, key in enumerate(actual):
            assert conditions[j], error_message(f"{key} FOR STATE {i}", expected[key], actual[key])


if __name__ == "__main__":
    test_quantum_dynamics_tully_1_real()
