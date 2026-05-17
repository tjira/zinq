from typing import Optional

from pydantic import BaseModel
from .surface_hopping import LandauZener

from ..potential import PotentialOptions


class LandauZenerOptions(BaseModel):
    def create(self, seed: int):
        return LandauZener(seed)


class SurfaceHoppingOptions(BaseModel):
    landau_zener: Optional[LandauZenerOptions] = None

    def create(self, seed: int):
        try: return next(opt.create(seed) for _, opt in self if opt is not None)
        except StopIteration: return None


class InitialConditionsOptions(BaseModel):
    position: list[float]
    momentum: list[float]
    gamma: list[float]
    state: int = 0


class LogIntervalOptions(BaseModel):
    iteration: int = 1
    trajectory: int = 1


class WriteOptions(BaseModel):
    kinetic_energy: Optional[str] = None
    momentum: Optional[str] = None
    population: Optional[str] = None
    position: Optional[str] = None
    potential_energy: Optional[str] = None
    total_energy: Optional[str] = None


class Options(BaseModel):
    initial_conditions: InitialConditionsOptions
    potential: PotentialOptions
    surface_hopping: Optional[SurfaceHoppingOptions] = None
    iterations: int
    trajectories: int = 1
    time_step: float
    mass: float = 1
    seed: int = 1
    log_interval: int = 1
    write: WriteOptions = WriteOptions()
