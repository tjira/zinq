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
        "trajectories": 100
    }

    result = run(options)

    expected = {
        "total_energy": 0.06594570188006,
        "kinetic_energy": 0.06614570217004,
        "potential_energy": -0.00020000028998,
        "position": [13.44416086839326],
        "momentum": [16.17365014916824],
        "population": [0.51000000000000, 0.49000000000000],
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


def test_classical_dynamics_harmonic_1d():
    options = {
        "initial_conditions": {
            "momentum": [0],
            "position": [1],
            "gamma": [2]
        },
        "potential": {
            "type" : "harmonic",
            "k": [1]
        },
        "log_interval": 200,
        "iterations": 1000,
        "time_step": 0.01,
        "trajectories": 100
    }

    result = run(options)

    expected = {
        "total_energy": 1.04931953645547,
        "kinetic_energy": 0.55008577542921,
        "potential_energy": 0.49923376102626,
        "position": [-0.76806238038594],
        "momentum": [0.58587360033957],
        "population": [1.00000000000000],
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


def test_classical_dynamics_harmonic_2d():
    options = {
        "initial_conditions": {
            "momentum": [0, 0],
            "position": [1, 1],
            "gamma": [2, 2]
        },
        "potential": {
            "type" : "harmonic",
            "k": [1, 1]
        },
        "log_interval": 200,
        "iterations": 1000,
        "time_step": 0.01,
        "trajectories": 100
    }

    result = run(options)

    expected = {
        "total_energy": 1.95121929902807,
        "kinetic_energy": 0.97203153303554,
        "potential_energy": 0.97918776599252,
        "position": [-0.69908887656263, -0.81874579090858],
        "momentum": [0.69180824756439, 0.50805324432117],
        "population": [1.00000000000000],
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


def test_classical_dynamics_harmonic_3d():
    options = {
        "initial_conditions": {
            "momentum": [0, 0, 0],
            "position": [1, 1, 1],
            "gamma": [2, 2, 2]
        },
        "potential": {
            "type" : "harmonic",
            "k": [1, 1, 1]
        },
        "log_interval": 200,
        "iterations": 1000,
        "time_step": 0.01,
        "trajectories": 100
    }

    result = run(options)

    expected = {
        "total_energy": 3.08488252209862,
        "kinetic_energy": 1.41235567464053,
        "potential_energy": 1.67252684745808,
        "position": [-0.90459347783846, -0.90122536304669, -0.60649409318960],
        "momentum": [0.51179517685186, 0.39142486146741, 0.60725988856574],
        "population": [1.00000000000000],
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


def test_classical_dynamics_fewest_switches_tully_1():
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
            "type": "fewest_switches"
        },
        "log_interval": 500,
        "mass": 2000,
        "iterations": 3000,
        "time_step": 1,
        "trajectories": 100
    }

    result = run(options)

    expected = {
        "total_energy": 0.06594330834986,
        "kinetic_energy": 0.06354330847411,
        "potential_energy": 0.00239999987575,
        "position": [13.20360646343850],
        "momentum": [15.86915258293628],
        "population": [0.38000000000000, 0.62000000000000],
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
