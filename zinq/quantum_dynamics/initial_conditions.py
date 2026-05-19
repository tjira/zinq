from dataclasses import dataclass

import numpy as np


@dataclass(frozen=True, kw_only=True)
class InitialConditions:
    pos: np.ndarray
    mom: np.ndarray
    gamma: np.ndarray
    state: int = 0
    adia: bool = False
