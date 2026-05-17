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
    def jump(self, ensemble: Ensemble, H: Hamiltonian, dt: float, time: float) -> None:
        pass

    def _propose_jumps(self, probs: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
        jumps = self._rng.random(probs.shape[0])[:, np.newaxis] < np.cumsum(probs, axis=1)
        jump_mask = np.any(jumps, axis=1)
        return jump_mask, np.argmax(jumps, axis=1)[jump_mask]

    def _apply_jump(self, ensemble: Ensemble, H: Hamiltonian, time: float, probs: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
        states, p = ensemble.states.copy(), ensemble.p.copy()
        jump_mask, new_states = self._propose_jumps(probs)

        if not np.any(jump_mask):
            return states, p

        V = H.pot.eval_a(list(ensemble.r.T), time)
        idx = np.flatnonzero(jump_mask)
        dV = V[new_states, idx] - V[ensemble.states[idx], idx]

        p_sq = (p[idx] ** 2).sum(axis=1)
        factor_sq = 1 - 2 * H.mass * dV / p_sq
        success = factor_sq > 0

        states[idx[success]] = new_states[success]
        p[idx[success]] *= np.sqrt(factor_sq[success, None])

        return states, p
