from typing import Literal

from pydantic import BaseModel, ConfigDict

from . import LandauZener


class LandauZenerOptions(BaseModel):
    model_config = ConfigDict(extra="forbid")

    type: Literal["landau_zener"]

    def create(self, seed: int):
        return LandauZener(seed)
