from dataclasses import dataclass, field
from typing import Literal

import numpy as np

from .potential import Potential


@dataclass(frozen=True, kw_only=True)
class Harmonic(Potential):
    name: Literal["harmonic"] = "harmonic"

    k: list[float] = field(default_factory=lambda: [1])

    @property
    def ndim(self) -> int:
        return len(self.k)

    @property
    def nstate(self) -> int:
        return 1

    def eval_d(self, r: list[np.ndarray], time: float = 0):
        V = sum((0.5 * self.k[i] * r[i]**2 for i in range(self.ndim)), start=np.zeros_like(r[0]))
        
        return V[..., np.newaxis, np.newaxis]
