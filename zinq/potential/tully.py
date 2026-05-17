from dataclasses import dataclass

from ..backend import np
from .potential import Potential


@dataclass(frozen=True, kw_only=True)
class TullyFirst(Potential):
    A: float = 0.01
    B: float = 1.6
    C: float = 0.005
    D: float = 1

    @property
    def ndim(self) -> int:
        return 1

    @property
    def nstate(self) -> int:
        return 2

    def eval_d(self, r: list[np.ndarray], time: float = 0.0):
        V00 = np.sign(r[0]) * self.A * (1 - np.exp(-self.B * np.abs(r[0])))
        V01 = self.C * np.exp(-self.D * r[0]**2)
        V11 = -V00

        return np.array([[V00, V01], [V01, V11]])
