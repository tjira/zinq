from dataclasses import dataclass
from typing import Literal

import numpy as np

from .potential import Potential


@dataclass(frozen=True, kw_only=True)
class TimeLinear(Potential):
    name: Literal["time_linear"] = "time_linear"

    coupling: float = 2
    slope: float = 10

    @property
    def ndim(self) -> int:
        return 1

    @property
    def nstate(self) -> int:
        return 2

    @property
    def is_td(self) -> bool:
        return True

    def eval_d(self, r: list[np.ndarray], time: float = 0.0) -> np.ndarray:
        V00 = np.full_like(r[0], self.slope * (time - self.slope))
        V01 = np.full_like(r[0], self.coupling)
        V11 = np.full_like(r[0], -V00)

        return np.stack([
            np.stack([V00, V01], axis=-1),
            np.stack([V01, V11], axis=-1)],
        axis=-2)
