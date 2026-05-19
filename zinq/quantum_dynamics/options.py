from typing import Annotated, Literal, Union

from pydantic import BaseModel, ConfigDict, Field

from ..potential import TullyFirst

AnyPotential = Annotated[Union[
    TullyFirst
], Field(discriminator="name")]


class GridConfig(BaseModel):
    model_config = ConfigDict(extra="forbid")

    limits: list[tuple[float, float]]
    npoint: int


class HamiltonianConfig(BaseModel):
    model_config = ConfigDict(extra="forbid")

    potential: AnyPotential
    mass: float


class ImaginaryConfig(BaseModel):
    model_config = ConfigDict(extra="forbid")

    nstate: int


class InitialConditionsConfig(BaseModel):
    model_config = ConfigDict(extra="forbid")

    position: list[float]
    momentum: list[float]
    gamma: list[float]
    state: int = 0
    adiabatic: bool = False


class Options(BaseModel):
    model_config = ConfigDict(extra="forbid")
    type: Literal["quantum_dynamics"]

    adiabatic: bool = False
    grid: GridConfig
    hamiltonian: HamiltonianConfig
    initial_conditions: InitialConditionsConfig
    imaginary: ImaginaryConfig | None = None
    iterations: int | None = None
    log_interval: int = 1
    time_step: float
