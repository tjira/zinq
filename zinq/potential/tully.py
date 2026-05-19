from dataclasses import dataclass
from typing import Literal

import numpy as np

from .potential import Potential


@dataclass(frozen=True, kw_only=True)
class TullyFirst(Potential):
    name: Literal["tully_1"] = "tully_1"

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

    def eval_d(self, r: list[np.ndarray], time: float = 0.0) -> np.ndarray:
        V00 = np.sign(r[0]) * self.A * (1 - np.exp(-self.B * np.abs(r[0])))
        V01 = self.C * np.exp(-self.D * r[0]**2)
        V11 = -V00

        return np.stack([
            np.stack([V00, V01], axis=-1),
            np.stack([V01, V11], axis=-1)],
        axis=-2)
