from ..backend import np


class Ensemble:
    r: np.ndarray
    p: np.ndarray
    c: np.ndarray
    states: np.ndarray
    mass: float

    def __init__(self, r: np.ndarray, p: np.ndarray, gamma: np.ndarray, mass: float, state: int, ntraj: int, nstate: int, seed: int = 1):
        assert mass > 0, "MASS MUST BE POSITIVE"
        for g in gamma:
            assert g >= 0, "GAMMA MUST BE NON-NEGATIVE"

        self.mass = mass
        self.r = np.zeros((ntraj, r.shape[0]))
        self.p = np.zeros((ntraj, p.shape[0]))
        self.states = np.full(ntraj, state, dtype=int)

        self.c = np.zeros((ntraj, nstate), dtype=complex)
        self.c[:, state] = 1

        rng = np.random.default_rng(seed)

        z_r = rng.standard_normal(size=(ntraj, r.shape[0]))
        z_p = rng.standard_normal(size=(ntraj, p.shape[0]))

        stdev_r = np.array([1 / np.sqrt(2 * g) if g > 0 else 0 for g in gamma])
        stdev_p = np.sqrt(gamma / 2)

        self.r = r + stdev_r * z_r
        self.p = p + stdev_p * z_p

    @property
    def ntraj(self) -> int:
        return self.r.shape[0]

    @property
    def ndim(self) -> int:
        return self.r.shape[1]

    def ke(self) -> float:
        return float(np.mean(np.sum(self.p**2, axis=1) / (2 * self.mass)))

    def pe(self, pot, time: float = 0) -> float:
        V = pot.eval_a(list(self.r.T), time)
        return float(np.mean(V[np.arange(self.ntraj), self.states]))

    def population(self, nstate: int) -> np.ndarray:
        return np.bincount(self.states, minlength=nstate) / self.ntraj

    def position(self) -> np.ndarray:
        return np.mean(self.r, axis=0)

    def momentum(self) -> np.ndarray:
        return np.mean(self.p, axis=0)
