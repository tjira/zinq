from dataclasses import dataclass

from ..backend import np


@dataclass(frozen=True)
class StateResult:
    kinetic_energy: float
    momentum: np.ndarray
    norm: float
    population: np.ndarray
    position: np.ndarray
    potential_energy: float
    total_energy: float

@dataclass(frozen=True)
class RunResult:
    states: list[StateResult]
