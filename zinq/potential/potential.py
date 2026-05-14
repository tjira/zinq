from abc import ABC, abstractmethod

import numpy as np

class Potential(ABC):
    @property
    @abstractmethod
    def ndim(self) -> int:
        pass

    @property
    @abstractmethod
    def nstate(self) -> int:
        pass

    @abstractmethod
    def evaluateDiabatic(self, r: list[np.ndarray], time: float = 0) -> np.ndarray:
        pass
