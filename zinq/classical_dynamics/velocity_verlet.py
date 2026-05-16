from typing import Optional

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
        return -self._numerical_gradient(position, states, time) / self.mass

    def _numerical_gradient(self, position: np.ndarray, states: np.ndarray, time: float) -> np.ndarray:
        ntraj, ndim, grad = *position.shape, np.zeros(position.shape)

        indices = np.arange(ntraj)

        for i in range(ndim):
            pos_plus = position.copy()
            pos_plus[:, i] += self.h
            
            pos_minus = position.copy()
            pos_minus[:, i] -= self.h
            
            r_plus = [pos_plus[:, j] for j in range(ndim)]
            r_minus = [pos_minus[:, j] for j in range(ndim)]
            
            V_plus = self.potential.eval_d(r_plus, time=time)
            V_minus = self.potential.eval_d(r_minus, time=time)
            
            grad[:, i] = (V_plus[states, states, indices] - V_minus[states, states, indices]) / (2 * self.h)
            
        return grad

    def step(self, ensemble: Ensemble, time: float) -> None:
        a_t = self._acceleration(ensemble.position, ensemble.states, time)
        ensemble.momentum += 0.5 * a_t * self.mass * self.dt
        ensemble.position += (ensemble.momentum / self.mass) * self.dt
        a_t_next = self._acceleration(ensemble.position, ensemble.states, time + self.dt)
        ensemble.momentum += 0.5 * a_t_next * self.mass * self.dt
