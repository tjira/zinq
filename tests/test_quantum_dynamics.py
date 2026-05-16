import json

import pytest

from zinq.quantum_dynamics.run import run


def error_message(observable: str, expected, actual) -> str:
    return f"INCORRECT {observable}: EXPECTED {expected}, GOT {actual}"


def test_quantum_dynamics_harmonic_1d_itp():
    options = {
        "grid": {
            "limits": [[-8, 8]],
            "npoint": 256
        },
        "initial_conditions": {
            "momentum": [0],
            "position": [1],
            "gamma": [2]
        },
        "potential": {
            "harmonic": {
                "k": [1]
            }
        },
        "log_interval": 200,
        "imaginary": {
            "nstate": 3
        },
        "iterations": 1000,
        "time_step": 0.01
    }

    result = run(options)

    expected_total_energies     = [0.50000000187134, 1.49999999864282, 2.49999999986185]
    expected_kinetic_energies   = [0.25000312532402, 0.75000937489134, 1.25001562583497]
    expected_potential_energies = [0.24999687654732, 0.74999062375148, 1.24998437402687]
    expected_norm_values        = [1.00000000000000, 1.00000000000000, 1.00000000000000]

    expected_position_values = [
        [ 0.00006053550970],
        [-0.00002269860870],
        [-0.00004994213660]
    ]
    expected_momentum_values = [
        [0.00000000000000],
        [0.00000000000000],
        [0.00000000000000]
    ]

    expected_populations = [
        [1.00000000000000],
        [1.00000000000000],
        [1.00000000000000]
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
        
        conditions = [actual[key] == pytest.approx(expected[key], abs=1e-12) for key in actual]

        for j, key in enumerate(actual):
            assert conditions[j], error_message(f"{key} FOR STATE {i}", expected[key], actual[key])


if __name__ == "__main__":
    test_quantum_dynamics_harmonic_1d_itp()
