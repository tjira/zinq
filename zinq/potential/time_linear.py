from dataclasses import dataclass

from pydantic import BaseModel

from ..backend import np
from .potential import Potential


class TimeLinearOptions(BaseModel):
    coupling: float = 2.0
    slope: float = 10.0

    def create(self) -> "TimeLinear":
        return TimeLinear(coupling=self.coupling, slope=self.slope)


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

    def eval_d(self, r: list[np.ndarray], time: float) -> np.ndarray:
        V00 = self.slope * (time - self.slope)
        V11 = -V00
        V01 = self.coupling

        return np.array([[V00, V01], [V01, V11]])
