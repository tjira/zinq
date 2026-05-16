from ..backend import np
from ..potential import Potential


class Ensemble:
    r: np.ndarray
    p: np.ndarray
    states: np.ndarray

    def __init__(self, pos: np.ndarray, mom: np.ndarray, gamma: np.ndarray, ntraj: int, state: int = 0, seed: int = 1):
        self.r = np.zeros((ntraj, pos.shape[0]))
        self.p = np.zeros((ntraj, mom.shape[0]))
        self.states = np.full(ntraj, state, dtype=int)

        stdev_pos, stdev_mom = 1 / np.sqrt(2 * gamma), np.sqrt(gamma / 2)

        rng = np.random.default_rng(seed)

        z_pos = rng.normal(size=(ntraj, pos.shape[0]))
        z_mom = rng.normal(size=(ntraj, mom.shape[0]))

        self.r = pos + stdev_pos * z_pos
        self.p = mom + stdev_mom * z_mom

    @property
    def ntraj(self) -> int:
        return self.r.shape[0]

    @property
    def ndim(self) -> int:
        return self.r.shape[1]

    def ke(self, mass: float) -> float:
        return float(np.mean(np.sum(self.p**2, axis=1) / (2 * mass)))

    def pe(self, potential: Potential, time: float = 0) -> float:
        v = potential.eval_a([self.r[:, j] for j in range(self.ndim)], time=time)
        return float(np.mean(v[self.states, np.arange(self.ntraj)]))

    def population(self, nstate: int) -> np.ndarray:
        return np.bincount(self.states, minlength=nstate) / self.ntraj

    def position(self) -> np.ndarray:
        return np.mean(self.r, axis=0)

    def momentum(self) -> np.ndarray:
        return np.mean(self.p, axis=0)
