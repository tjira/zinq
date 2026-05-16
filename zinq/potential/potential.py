from abc import ABC, abstractmethod

from ..backend import np


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
    def eval_d(self, r: list[np.ndarray], time: float = 0) -> np.ndarray:
        pass

    def eval_a(self, r: list[np.ndarray], time: float = 0) -> np.ndarray:
        v_d = self.eval_d(r, time)

        if v_d.ndim == 2: return np.linalg.eigvalsh(v_d)

        return np.linalg.eigvalsh(v_d.transpose(2, 0, 1)).transpose()
