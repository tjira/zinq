from dataclasses import dataclass

import numpy as np


@dataclass(frozen=True, kw_only=True)
class State:
    total_energy: float
    kinetic_energy: float
    potential_energy: float
    position: np.ndarray
    momentum: np.ndarray
    norm: float
    population: np.ndarray


@dataclass(frozen=True, kw_only=True)
class Results:
    states: list[State]
