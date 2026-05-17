from typing import Literal

from pydantic import BaseModel

from . import LandauZener


class LandauZenerOptions(BaseModel):
    type: Literal["landau_zener"]

    def create(self, seed: int):
        return LandauZener(seed)
