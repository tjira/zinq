from typing import Literal, Optional, Union

from pydantic import BaseModel, Field

from ..potential import HarmonicOptions, TimeLinearOptions, TullyFirstOptions
from .surface_hopping import LandauZener, LandauZenerOptions


class InitialConditionsOptions(BaseModel):
    gamma: list[float]
    momentum: list[float]
    position: list[float]
    state: int = 0


class WriteOptions(BaseModel):
    kinetic_energy: Optional[str] = None
    momentum: Optional[str] = None
    population: Optional[str] = None
    position: Optional[str] = None
    potential_energy: Optional[str] = None
    total_energy: Optional[str] = None


class Options(BaseModel):
    initial_conditions: InitialConditionsOptions
    iterations: int
    log_interval: int = 1
    mass: float = 1
    seed: int = 1
    time_step: float
    trajectories: int = 1
    write: WriteOptions = WriteOptions()
    potential: Union[
        HarmonicOptions,
        TimeLinearOptions,
        TullyFirstOptions,
    ] = Field(discriminator="type")
    surface_hopping: Optional[Union[
        LandauZenerOptions,
    ]] = Field(default=None, discriminator="type")
