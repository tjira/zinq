from dataclasses import dataclass

from ..backend import np


@dataclass(frozen=True)
class RunResult:
    kinetic_energy: float
    momentum: np.ndarray
    population: np.ndarray
    position: np.ndarray
    potential_energy: float
    total_energy: float
