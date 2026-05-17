from typing import Literal, Optional, Union

from pydantic import BaseModel, ConfigDict, Field

from ..potential import HarmonicOptions, TimeLinearOptions, TullyFirstOptions
from .surface_hopping import FewestSwitchesOptions, LandauZenerOptions


class InitialConditionsOptions(BaseModel):
    model_config = ConfigDict(extra="forbid")

    gamma: list[float]
    momentum: list[float]
    position: list[float]
    state: int = 0


class WriteOptions(BaseModel):
    model_config = ConfigDict(extra="forbid")

    kinetic_energy: Optional[str] = None
    momentum: Optional[str] = None
    population: Optional[str] = None
    position: Optional[str] = None
    potential_energy: Optional[str] = None
    total_energy: Optional[str] = None


class Options(BaseModel):
    model_config = ConfigDict(extra="forbid")

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
        FewestSwitchesOptions,
    ]] = Field(default=None, discriminator="type")
