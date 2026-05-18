from typing import Optional
from abc import ABC, abstractmethod

from ...backend import np
from ..ensemble import Ensemble
from ...potential import Potential


class SurfaceHopping(ABC):
    _rng: np.random.Generator

    def __init__(self, seed: int) -> None:
        self._rng = np.random.default_rng(seed)

    @abstractmethod
    def jump(self, ensemble: Ensemble, pot: Potential, dt: float, time: float) -> None:
        pass

    def _propose_jumps(self, probs: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
        jumps = self._rng.random(probs.shape[0])[:, np.newaxis] < np.cumsum(probs, axis=1)
        jump_mask = np.any(jumps, axis=1)
        return jump_mask, np.argmax(jumps, axis=1)[jump_mask]

    def _apply_jump(self, ensemble: Ensemble, V_a: np.ndarray, probs: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
        states, p = ensemble.states.copy(), ensemble.p.copy()
        jump_mask, new_states = self._propose_jumps(probs)

        if not np.any(jump_mask):
            return states, p

        idx = np.flatnonzero(jump_mask)
        dV = V_a[idx, new_states] - V_a[idx, ensemble.states[idx]]

        p_sq = (p[idx] ** 2).sum(axis=1)

        mask_nonzero = p_sq > 1e-14
        factor_sq = np.zeros_like(p_sq)
        factor_sq[mask_nonzero] = 1 - 2 * ensemble.mass * dV[mask_nonzero] / p_sq[mask_nonzero]

        success = (factor_sq > 0) & mask_nonzero

        states[idx[success]] = new_states[success]
        p[idx[success]] *= np.sqrt(factor_sq[success, None])

        return states, p
