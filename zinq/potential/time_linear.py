from dataclasses import dataclass

from ..backend import np
from .potential import Potential


@dataclass(frozen=True, kw_only=True)
class TimeLinear(Potential):
    coupling: float = 2.0
    slope: float = 10.0

    @property
    def ndim(self) -> int:
        return 1

    @property
    def nstate(self) -> int:
        return 2

    @property
    def is_time_dependent(self) -> bool:
        return True

    def eval_d(self, r: list[np.ndarray], time: float = 0.0) -> np.ndarray:
        V00 = np.full_like(r[0], self.slope * (time - self.slope))
        V01 = np.full_like(r[0], self.coupling)
        V11 = np.full_like(r[0], -V00)

        return np.moveaxis(np.array([[V00, V01], [V01, V11]]), [0, 1], [-2, -1])
