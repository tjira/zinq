from dataclasses import dataclass

import numpy as np


@dataclass(frozen=True, kw_only=True)
class Result:
    total_energy: np.ndarray
    kinetic_energy: np.ndarray
    potential_energy: np.ndarray
    position: np.ndarray
    momentum: np.ndarray
    norm: np.ndarray
    population: np.ndarray
