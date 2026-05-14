from dataclasses import dataclass

import numpy as np

from pydantic import BaseModel

from .potential import Potential

class HarmonicOptions(BaseModel):
    k: list[float] = [1]

    def create(self) -> "Harmonic":
        return Harmonic(k=np.array(self.k))

@dataclass(frozen=True, kw_only=True)
class Harmonic(Potential):
    k: np.ndarray

    @property
    def ndim(self) -> int:
        return self.k.shape[0]

    @property
    def nstate(self) -> int:
        return 1

    def evaluateDiabatic(self, r: list[np.ndarray], time: float = 0):
        assert len(r) == self.ndim
        
        V = sum(0.5 * self.k[i] * r[i]**2 for i in range(self.ndim))
        
        return V[np.newaxis, np.newaxis, ...]
