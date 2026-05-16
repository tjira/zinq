from dataclasses import dataclass

from ..backend import np


@dataclass(frozen=True)
class RunResult:
    population: np.ndarray
    kinetic_energy: float
    potential_energy: float
    total_energy: float
    position: np.ndarray
    momentum: np.ndarray
