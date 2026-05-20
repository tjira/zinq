from abc import ABC, abstractmethod
from typing import Annotated, Union

import numpy as np
from pydantic import Field


class Potential(ABC):
    @property
    @abstractmethod
    def ndim(self) -> int:
        pass

    @property
    @abstractmethod
    def nstate(self) -> int:
        pass

    @property
    def is_td(self) -> bool:
        return False

    @abstractmethod
    def eval_d(self, r: list[np.ndarray], time: float) -> np.ndarray:
        pass

    def eval_a(self, r: list[np.ndarray], time: float) -> np.ndarray:
        return np.linalg.eigvalsh(self.eval_d(r, time))


from .harmonic import Harmonic
from .tully import TullyFirst

AnyPotential = Annotated[Union[
    Harmonic,
    TullyFirst,
], Field(discriminator="name")]
