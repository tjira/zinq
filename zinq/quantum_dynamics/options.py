from typing import Optional

from pydantic import BaseModel

from ..potential import HarmonicOptions, Potential, TullyFirstOptions


class GridOptions(BaseModel):
    limits: list[list[float]]
    npoint: int

class ImaginaryOptions(BaseModel):
    nstate: int = 1

class InitialConditionsOptions(BaseModel):
    position: list[float]
    momentum: list[float]
    gamma: list[float]
    state: int = 0

class PotentialOptions(BaseModel):
    harmonic: Optional[HarmonicOptions] = None
    tully_1: Optional[TullyFirstOptions] = None

    def create(self) -> Potential:
        try: return next(pot.create() for _, pot in self if pot is not None)
        except StopIteration: raise ValueError("NO POTENTIAL SPECIFIED")

class Options(BaseModel):
    grid: GridOptions
    initial_conditions: InitialConditionsOptions
    potential: PotentialOptions
    iterations: int
    time_step: float
    mass: float = 1
    log_interval: int = 1
    imaginary: Optional[ImaginaryOptions] = None
