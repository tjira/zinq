from pydantic import BaseModel

from .surface_hopping import SurfaceHopping
from ...backend import np
from ...potential import Potential
from ..ensemble import Ensemble


class LandauZenerOptions(BaseModel):
    def create(self):
        return LandauZener()


class LandauZener(SurfaceHopping):
    def jump(self, ensemble: Ensemble, potential: Potential, dt: float, time: float) -> None:
        pass
