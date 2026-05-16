from typing import Optional

from pydantic import BaseModel

from .harmonic import HarmonicOptions
from .potential import Potential
from .time_linear import TimeLinearOptions
from .tully import TullyFirstOptions


class PotentialOptions(BaseModel):
    harmonic: Optional[HarmonicOptions] = None
    time_linear: Optional[TimeLinearOptions] = None
    tully_1: Optional[TullyFirstOptions] = None

    def create(self) -> Potential:
        try: return next(pot.create() for _, pot in self if pot is not None)
        except StopIteration: raise ValueError("NO POTENTIAL SPECIFIED")
