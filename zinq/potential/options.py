from typing import Literal

from pydantic import BaseModel

from ..backend import np
from .harmonic import Harmonic
from .time_linear import TimeLinear
from .tully import TullyFirst


class HarmonicOptions(BaseModel):
    type: Literal["harmonic"]

    k: list[float] = [1]

    def create(self):
        return Harmonic(k=np.array(self.k))


class TimeLinearOptions(BaseModel):
    type: Literal["time_linear"]

    coupling: float = 2
    slope: float = 10

    def create(self):
        return TimeLinear(coupling=self.coupling, slope=self.slope)


class TullyFirstOptions(BaseModel):
    type: Literal["tully_1"]

    A: float = 0.01
    B: float = 1.6
    C: float = 0.005
    D: float = 1.0

    def create(self):
        return TullyFirst(A=self.A, B=self.B, C=self.C, D=self.D)
