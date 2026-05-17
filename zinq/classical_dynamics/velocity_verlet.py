from typing import Callable, Optional

from ..backend import np
from .ensemble import Ensemble
from .hamiltonian import Hamiltonian


class VelocityVerlet:
    H: Hamiltonian
    dt: float
    h: float

    def __init__(self, H: Hamiltonian, dt: float, h: float = 1e-8):
        self.H, self.dt, self.h = H, dt, h

    def _acceleration(self, position: np.ndarray, states: np.ndarray, time: float) -> np.ndarray:
        grad = self.H.pot.grad_a(list(position.T), time, self.h)

        return -grad[:, states, np.arange(position.shape[0])].T / self.H.mass

    def step(self, ensemble: Ensemble, time: float, jump_fn: Optional[Callable] = None) -> None:
        a_t = self._acceleration(ensemble.r, ensemble.states, time)
        ensemble.p += 0.5 * a_t * self.H.mass * self.dt
        ensemble.r += (ensemble.p / self.H.mass) * self.dt

        if jump_fn:
            jump_fn(ensemble, self.H, self.dt, time + self.dt)

        a_t_next = self._acceleration(ensemble.r, ensemble.states, time + self.dt)
        ensemble.p += 0.5 * a_t_next * self.H.mass * self.dt
