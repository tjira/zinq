import pytest

from zinq.classical_dynamics.run import run


def error_message(observable: str, expected, actual) -> str:
    return f"INCORRECT {observable}: EXPECTED {expected}, GOT {actual}"

def test_classical_dynamics_landau_zener_tully_1():
    options = {
        "initial_conditions": {
            "momentum": [15],
            "gamma": [2],
            "position": [-10],
            "state": 1
        },
        "potential": {
            "type": "tully_1"
        },
        "surface_hopping": {
            "type": "landau_zener"
        },
        "log_interval": 500,
        "mass": 2000,
        "iterations": 3000,
        "time_step": 1,
        "trajectories": 10000
    }

    result = run(options)

    expected = {
        "total_energy": 0.06640837161028,
        "kinetic_energy": 0.06630837178705,
        "potential_energy": 0.00009999982323,
        "position": [13.53720742261733],
        "momentum": [16.21021096251392],
        "population": [0.49500000000000, 0.50500000000000],
    }

    actual = {
        "total_energy": result.total_energy,
        "kinetic_energy": result.kinetic_energy,
        "potential_energy": result.potential_energy,
        "position": result.position.tolist(),
        "momentum": result.momentum.tolist(),
        "population": result.population.tolist(),
    }

    conditions = [actual[key] == pytest.approx(expected[key], abs=1e-08) for key in actual]

    for j, key in enumerate(actual):
        assert conditions[j], error_message(f"{key}", expected[key], actual[key])
