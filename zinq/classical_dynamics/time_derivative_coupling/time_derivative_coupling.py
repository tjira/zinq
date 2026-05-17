from abc import ABC, abstractmethod

from ...backend import np


class TimeDerivativeCoupling(ABC):
    @abstractmethod
    def evaluate(self, U_old: np.ndarray, U_new: np.ndarray, dt: float) -> np.ndarray:
        pass
