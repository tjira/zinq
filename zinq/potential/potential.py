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
        return np.linalg.eigvalsh(self.eval_d(r, time))

    def grad_a(self, r: list[np.ndarray], time: float = 0.0, step: float = 1e-8) -> np.ndarray:
        grad = []

        for i in range(len(r)):
            r_plus = r.copy()
            r_minus = r.copy()

            r_plus[i] = r_plus[i] + step
            r_minus[i] = r_minus[i] - step

            V_plus = self.eval_a(r_plus, time)
            V_minus = self.eval_a(r_minus, time)

            grad.append((V_plus - V_minus) / (2 * step))

        return np.moveaxis(np.array(grad), 0, -2)
