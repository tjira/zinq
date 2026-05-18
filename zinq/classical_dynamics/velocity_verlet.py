from typing import Callable, Optional

from ..backend import np
from .ensemble import Ensemble
from ..potential import Potential


class VelocityVerlet:
    pot: Potential
    dt: float
    h: float

    def __init__(self, pot: Potential, dt: float, h: float = 1e-8):
        self.pot, self.dt, self.h = pot, dt, h

    def _acceleration(self, position: np.ndarray, states: np.ndarray, mass: float, time: float) -> np.ndarray:
        grad = self.pot.grad_a(list(position.T), time, self.h)
        return -grad[np.arange(position.shape[0]), :, states] / mass

    def step(self, ensemble: Ensemble, time: float, jump_fn: Optional[Callable] = None) -> None:
        a_t = self._acceleration(ensemble.r, ensemble.states, ensemble.mass, time)
        ensemble.p += 0.5 * a_t * ensemble.mass * self.dt
        ensemble.r += (ensemble.p / ensemble.mass) * self.dt

        if jump_fn: jump_fn(ensemble, self.pot, self.dt, time + self.dt)

        a_t_next = self._acceleration(ensemble.r, ensemble.states, ensemble.mass, time + self.dt)
        ensemble.p += 0.5 * a_t_next * ensemble.mass * self.dt
