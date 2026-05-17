from typing import Optional
from abc import ABC, abstractmethod

from ...backend import np
from ..ensemble import Ensemble
from ..hamiltonian import Hamiltonian


class SurfaceHopping(ABC):
    _rng: np.random.Generator

    def __init__(self, seed: int) -> None:
        self._rng = np.random.default_rng(seed)

    @abstractmethod
    def probabilities(self, ensemble: Ensemble, H: Hamiltonian, dt: float, time: float) -> Optional[np.ndarray]:
        pass

    def _propose_jumps(self, probs: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
        jumps = self._rng.random(probs.shape[0])[:, np.newaxis] < np.cumsum(probs, axis=1)
        
        return np.any(jumps, axis=1), np.argmax(jumps, axis=1)[np.any(jumps, axis=1)]

    def jump(self, ensemble: Ensemble, H: Hamiltonian, dt: float, time: float):
        if (probs := self.probabilities(ensemble, H, dt, time)) is None: return

        jump_mask, new_states = self._propose_jumps(probs)

        if not np.any(jump_mask): return

        V, idx = H.pot.eval_a(list(ensemble.r.T), time), np.flatnonzero(jump_mask)

        dV = V[new_states, idx := np.flatnonzero(jump_mask)] - V[ensemble.states[idx], idx]

        p_sq = (ensemble.p[idx] ** 2).sum(axis=1)
        factor_sq = 1 - 2 * H.mass * dV / p_sq

        success = factor_sq > 0

        ensemble.states[idx[success]] = new_states[success]
        ensemble.p[idx[success]] *= np.sqrt(factor_sq[success, None])
