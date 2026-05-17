from abc import ABC, abstractmethod

from ...backend import np
from ..ensemble import Ensemble
from ..hamiltonian import Hamiltonian


class SurfaceHopping(ABC):
    @abstractmethod
    def jump(self, ensemble: Ensemble, H: Hamiltonian, dt: float, time: float) -> None:
        pass
