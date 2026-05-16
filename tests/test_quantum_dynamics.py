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
    expected_potential_energies = [0.24999687654732, 0.74999062375147, 1.24998437402687]
    expected_norm_values        = [1.00000000000000, 1.00000000000000, 1.00000000000000]

    expected_position_values = [
        [ 0.00006053550965],
        [-0.00002269860869],
        [-0.00004994213663]
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
        
        conditions = [actual[key] == pytest.approx(expected[key], abs=1e-08) for key in actual]

        for j, key in enumerate(actual):
            assert conditions[j], error_message(f"{key} FOR STATE {i}", expected[key], actual[key])


def test_quantum_dynamics_harmonic_2d_itp():
    options = {
        "grid": {
            "limits": [[-8, 8], [-8, 8]],
            "npoint": 256
        },
        "initial_conditions": {
            "momentum": [0, 0],
            "position": [1, 1],
            "gamma": [2, 2]
        },
        "potential": {
            "harmonic": {
                "k": [1, 1]
            }
        },
        "log_interval": 200,
        "imaginary": {
            "nstate": 2
        },
        "iterations": 1000,
        "time_step": 0.01
    }

    result = run(options)

    expected_total_energies     = [1.00000000374269, 1.99999999776581]
    expected_kinetic_energies   = [0.50000625064804, 1.00001249941377]
    expected_potential_energies = [0.49999375309464, 0.99998749835204]
    expected_norm_values        = [1.00000000000000, 1.00000000000000]

    expected_position_values = [
        [ 0.00006053550965,  0.00006053550965],
        [-0.00001134791748, -0.00001134791748]
    ]
    expected_momentum_values = [
        [0.00000000000000, 0.00000000000000],
        [0.00000000000000, 0.00000000000000]
    ]

    expected_populations = [
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
        
        conditions = [actual[key] == pytest.approx(expected[key], abs=1e-08) for key in actual]

        for j, key in enumerate(actual):
            assert conditions[j], error_message(f"{key} FOR STATE {i}", expected[key], actual[key])


def test_quantum_dynamics_tully_1_real():
    options = {
        "grid": {
            "limits": [[-24, 24]],
            "npoint": 2048
        },
        "initial_conditions": {
            "momentum": [15],
            "position": [-10],
            "gamma": [2],
            "state": 1
        },
        "potential": {
            "tully_1": {}
        },
        "log_interval": 500,
        "mass": 2000,
        "iterations": 3000,
        "time_step": 1
    }

    result = run(options)

    expected_total_energies     = [0.06649999845028]
    expected_kinetic_energies   = [0.06470766411999]
    expected_potential_energies = [0.00179233433028]
    expected_norm_values        = [1.00000000000074]

    expected_position_values = [
        [13.41125194061655]
    ]

    expected_momentum_values = [
        [16.01020083373733]
    ]

    expected_populations = [
        [0.58961671754100, 0.41038328245899]
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
    test_quantum_dynamics_harmonic_1d_itp()
    test_quantum_dynamics_harmonic_2d_itp()
    test_quantum_dynamics_tully_1_real()
