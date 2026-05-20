from dataclasses import dataclass

import numpy as np


@dataclass(frozen=True, kw_only=True)
class Observables:
    norm: float
    pop: np.ndarray
    pos: np.ndarray
    mom: np.ndarray
    pe: float
    ke: float
    e: float
