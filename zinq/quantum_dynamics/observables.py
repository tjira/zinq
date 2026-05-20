from dataclasses import dataclass

import numpy as np


@dataclass(frozen=True, kw_only=True)
class Observables:
    acf: complex | None = None
    e: float | None = None
    ke: float | None = None
    mom: np.ndarray | None = None
    norm: float | None = None
    pe: float | None = None
    pop: np.ndarray | None = None
    pos: np.ndarray | None = None
    wfn: np.ndarray | None = None
