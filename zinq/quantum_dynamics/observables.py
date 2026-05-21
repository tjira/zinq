from dataclasses import dataclass

import numpy as np


@dataclass(frozen=True, kw_only=True)
class Observables:
    autocorrelation: complex | None = None
    total_energy: float | None = None
    kinetic_energy: float | None = None
    momentum: np.ndarray | None = None
    norm: float | None = None
    potential_energy: float | None = None
    population: np.ndarray | None = None
    position: np.ndarray | None = None
    wavefunction: np.ndarray | None = None
