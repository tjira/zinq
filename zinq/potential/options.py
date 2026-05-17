from typing import Optional

from pydantic import BaseModel
from .harmonic import Harmonic
from .time_linear import TimeLinear
from ..backend import np
from .tully import TullyFirst

from .potential import Potential


class HarmonicOptions(BaseModel):
    k: list[float] = [1]

    def create(self):
        return Harmonic(k=np.array(self.k))


class TimeLinearOptions(BaseModel):
    coupling: float = 2.0
    slope: float = 10.0

    def create(self):
        return TimeLinear(coupling=self.coupling, slope=self.slope)


class TullyFirstOptions(BaseModel):
    A: float = 0.01
    B: float = 1.6
    C: float = 0.005
    D: float = 1.0

    def create(self):
        return TullyFirst(A=self.A, B=self.B, C=self.C, D=self.D)


class PotentialOptions(BaseModel):
    harmonic: Optional[HarmonicOptions] = None
    time_linear: Optional[TimeLinearOptions] = None
    tully_1: Optional[TullyFirstOptions] = None

    def create(self) -> Potential:
        try: return next(pot.create() for _, pot in self if pot is not None)
        except StopIteration: raise ValueError("NO POTENTIAL SPECIFIED")
