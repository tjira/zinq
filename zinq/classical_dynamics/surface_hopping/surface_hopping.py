from abc import ABC, abstractmethod

from ...backend import np
from ...potential import Potential
from ..ensemble import Ensemble


class SurfaceHopping(ABC):
    @abstractmethod
    def jump(self, ensemble: Ensemble, pot: Potential, dt: float, time: float) -> None:
        pass
