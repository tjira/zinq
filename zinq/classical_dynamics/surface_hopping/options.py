from typing import Literal, Union

from pydantic import BaseModel, ConfigDict, Field

from ..time_derivative_coupling import HammesSchifferTullyOptions
from .fewest_switches import FewestSwitches
from .landau_zener import LandauZener


class LandauZenerOptions(BaseModel):
    model_config = ConfigDict(extra="forbid")

    type: Literal["landau_zener"]

    def create(self, seed: int):
        return LandauZener(seed)


class FewestSwitchesOptions(BaseModel):
    model_config = ConfigDict(extra="forbid")

    type: Literal["fewest_switches"]

    substeps: int = 10
    tdc: Union[
        HammesSchifferTullyOptions,
    ] = Field(HammesSchifferTullyOptions(type="hammes_schiffer_tully"), discriminator="type")

    def create(self, seed: int):
        return FewestSwitches(seed, self.tdc.create(), self.substeps)
