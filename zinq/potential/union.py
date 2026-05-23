"""Module defining the AnyPotential union type for concrete potentials."""

from typing import Annotated

from pydantic import Field

from .harmonic import Harmonic
from .time_linear import TimeLinear
from .tully import TullyFirst

AnyPotential = Annotated[Harmonic | TimeLinear | TullyFirst, Field(discriminator="name")]
