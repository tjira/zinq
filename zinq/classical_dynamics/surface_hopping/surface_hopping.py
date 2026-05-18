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

    def _eval_hops(self, p_proposed: np.ndarray, mass: float, dV: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
        p_sq = np.sum(p_proposed ** 2, axis=1)

        mask_nonzero = p_sq > 1e-14
        factor_sq = np.zeros_like(p_sq)
        factor_sq[mask_nonzero] = 1 - 2 * mass * dV[mask_nonzero] / p_sq[mask_nonzero]

        success = (factor_sq > 0) & mask_nonzero

        return success, factor_sq

    def _propose_jumps(self, probs: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
        jumps = self._rng.random(probs.shape[0])[:, np.newaxis] < np.cumsum(probs, axis=1)
        return np.any(jumps, axis=1), np.argmax(jumps, axis=1)

    def _apply_jump(self, ensemble: Ensemble, V_a: np.ndarray, jump_mask: np.ndarray, target_states: np.ndarray) -> np.ndarray:
        if not np.any(jump_mask):
            return np.zeros(ensemble.ntraj, dtype=bool)

        idx = np.flatnonzero(jump_mask)
        dV = V_a[idx, target_states[idx]] - V_a[idx, ensemble.states[idx]]

        success, factor_sq = self._eval_hops(ensemble.p[idx], ensemble.mass, dV)

        ensemble.states[idx[success]] = target_states[idx[success]]
        ensemble.p[idx[success]] *= np.sqrt(factor_sq[success, None])

        return np.bincount(idx[success], minlength=ensemble.ntraj).astype(bool)
