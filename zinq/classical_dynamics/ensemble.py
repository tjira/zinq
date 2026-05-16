from ..backend import np


class Ensemble:
    position: np.ndarray
    momentum: np.ndarray
    states: np.ndarray

    def __init__(self, position: np.ndarray, momentum: np.ndarray, gamma: np.ndarray, ntraj: int, state: int = 0):
        self.position = np.zeros((ntraj, position.shape[0]))
        self.momentum = np.zeros((ntraj, momentum.shape[0]))
        self.states = np.full(ntraj, state, dtype=int)

        stdev_pos, stdev_mom = 1 / np.sqrt(2 * gamma), np.sqrt(gamma / 2)

        z_pos = np.random.normal(size=(ntraj, position.shape[0]))
        z_mom = np.random.normal(size=(ntraj, momentum.shape[0]))

        self.position = position + stdev_pos * z_pos
        self.momentum = momentum + stdev_mom * z_mom

    @property
    def ntraj(self) -> int:
        return self.position.shape[0]

    @property
    def ndim(self) -> int:
        return self.position.shape[1]

    def ke(self) -> np.ndarray:
        return np.mean(np.sum(self.momentum**2, axis=1) / 2)

    def pe(self, potential) -> float:
        V = potential.eval_d([self.position[:, j] for j in range(self.ndim)])
        return np.mean(V[self.states, self.states, np.arange(self.ntraj)])
