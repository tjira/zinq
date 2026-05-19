from typing import Annotated, Literal, Union

from pydantic import BaseModel, ConfigDict, Field

from ..potential import TullyFirst
from .grid import Grid
from .initial_conditions import InitialConditions

AnyPotential = Annotated[Union[
    TullyFirst
], Field(discriminator="name")]


class Options(BaseModel):
    model_config = ConfigDict(extra="forbid")
    type: Literal["quantum_dynamics"]

    adiabatic: bool = False
    grid: Grid.Config
    initial_conditions: InitialConditions.Config
    iterations: int | None = None
    log_interval: int = 1
    mass: float
    potential: AnyPotential
    time_step: float
