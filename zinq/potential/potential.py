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

    @property
    def is_time_dependent(self) -> bool:
        return False

    @abstractmethod
    def eval_d(self, r: list[np.ndarray], time: float) -> np.ndarray:
        pass

    def eval_a(self, r: list[np.ndarray], time: float) -> np.ndarray:
        return np.linalg.eigvalsh(self.eval_d(r, time))
