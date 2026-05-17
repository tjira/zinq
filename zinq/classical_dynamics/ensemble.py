from ..backend import np
from .hamiltonian import Hamiltonian


class Ensemble:
    r: np.ndarray
    p: np.ndarray
    states: np.ndarray

    def __init__(self, r: np.ndarray, p: np.ndarray, gamma: np.ndarray, ntraj: int, state: int = 0, seed: int = 1):
        self.r = np.zeros((ntraj, r.shape[0]))
        self.p = np.zeros((ntraj, p.shape[0]))

        self.states = np.full(ntraj, state, dtype=int)

        stdev_r, stdev_p = 1 / np.sqrt(2 * gamma), np.sqrt(gamma / 2)

        rng = np.random.default_rng(seed)

        z_r = rng.normal(size=(ntraj, r.shape[0]))
        z_p = rng.normal(size=(ntraj, p.shape[0]))

        self.r = r + stdev_r * z_r
        self.p = p + stdev_p * z_p

    @property
    def ntraj(self) -> int:
        return self.r.shape[0]

    @property
    def ndim(self) -> int:
        return self.r.shape[1]

    def ke(self, H: Hamiltonian) -> float:
        return float(np.mean(np.sum(self.p**2, axis=1) / (2 * H.mass)))

    def pe(self, H: Hamiltonian, time: float = 0) -> float:
        V = H.pot.eval_a([self.r[:, j] for j in range(self.ndim)], time=time)
        return float(np.mean(V[self.states, np.arange(self.ntraj)]))

    def population(self, nstate: int) -> np.ndarray:
        return np.bincount(self.states, minlength=nstate) / self.ntraj

    def position(self) -> np.ndarray:
        return np.mean(self.r, axis=0)

    def momentum(self) -> np.ndarray:
        return np.mean(self.p, axis=0)
