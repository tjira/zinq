from dataclasses import dataclass
from ..backend import np

@dataclass(frozen=True)
class StateResult:
    population: np.ndarray
    kinetic_energy: float
    potential_energy: float
    total_energy: float
    position: np.ndarray
    momentum: np.ndarray
    norm: float

@dataclass(frozen=True)
class RunResult:
    states: list[StateResult]
