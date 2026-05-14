from pydantic import BaseModel
from typing import Optional

from ..potential import HarmonicOptions, TullyFirstOptions

class Grid(BaseModel):
    limits: list[list[float]]
    npoint: int

class InitialConditions(BaseModel):
    position: list[float]
    momentum: list[float]
    gamma: list[float]
    state: int = 0

class PotentialOptions(BaseModel):
    harmonic: Optional[HarmonicOptions] = None
    tully_1: Optional[TullyFirstOptions] = None

class Options(BaseModel):
    grid: Grid
    initial_conditions: InitialConditions
    potential: PotentialOptions
    mass: float
    iterations: int
    time_step: float
    log_interval: int = 1
