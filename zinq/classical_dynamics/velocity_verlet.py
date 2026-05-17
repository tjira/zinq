from typing import Callable, Optional

from ..backend import np
from ..potential import Potential
from .ensemble import Ensemble


class VelocityVerlet:
    potential: Potential
    mass: float
    dt: float
    h: float

    def __init__(self, potential: Potential, mass: float, dt: float, h: float = 1e-8):
        self.potential, self.mass, self.dt, self.h = potential, mass, dt, h

    def _acceleration(self, position: np.ndarray, states: np.ndarray, time: float) -> np.ndarray:
        grad = self.potential.grad_a(list(position.T), time, self.h)

        return -grad[:, states, np.arange(position.shape[0])].T / self.mass

    def step(self, ensemble: Ensemble, time: float, jump_fn: Optional[Callable] = None) -> None:
        a_t = self._acceleration(ensemble.r, ensemble.states, time)
        ensemble.p += 0.5 * a_t * self.mass * self.dt
        ensemble.r += (ensemble.p / self.mass) * self.dt

        if jump_fn:
            jump_fn(ensemble, self.potential, self.dt, time + self.dt)

        a_t_next = self._acceleration(ensemble.r, ensemble.states, time + self.dt)
        ensemble.p += 0.5 * a_t_next * self.mass * self.dt
