from dataclasses import dataclass

import numpy as np
from pydantic import BaseModel

from .potential import Potential


class TullyFirstOptions(BaseModel):
    A: float = 0.01
    B: float = 1.6
    C: float = 0.005
    D: float = 1.0

    def create(self):
        return TullyFirst(A=self.A, B=self.B, C=self.C, D=self.D)


@dataclass(frozen=True, kw_only=True)
class TullyFirst(Potential):
    ndim: int = 1
    nstate: int = 2

    A: float = 0.01
    B: float = 1.6
    C: float = 0.005
    D: float = 1

    def evaluateDiabatic(self, r: list[np.ndarray], time: float = 0.0):
        V00 = np.sign(r[0]) * self.A * (1 - np.exp(-self.B * np.abs(r[0])))
        V01 = self.C * np.exp(-self.D * r[0]**2)
        V22 = np.sign(r[0]) * self.A * (np.exp(-self.B * np.abs(r[0])) - 1)

        return np.array([[V00, V01], [V01, V22]])
