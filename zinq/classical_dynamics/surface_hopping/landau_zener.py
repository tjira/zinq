from typing import Optional
from collections import deque

from ...backend import np
from ..ensemble import Ensemble
from ..hamiltonian import Hamiltonian
from .surface_hopping import SurfaceHopping


class LandauZener(SurfaceHopping):
    _V_history: deque

    def __init__(self, seed: int) -> None:
        super().__init__(seed); self._V_history = deque(maxlen=3)

    def _calc_gaps(self, states: np.ndarray, state: int, dt: float) -> tuple[np.ndarray, ...]:
        indices = np.arange(states.shape[0])

        z0 = np.abs(self._V_history[2][:, state] - self._V_history[2][indices, states])
        z1 = np.abs(self._V_history[1][:, state] - self._V_history[1][indices, states])
        z2 = np.abs(self._V_history[0][:, state] - self._V_history[0][indices, states])

        dz0, dz1, ddz1 = (z0 - z1) / dt, (z1 - z2) / dt, (z0 - 2 * z1 + z2) / (dt * dt)

        return z0, z1, z2, dz0, dz1, ddz1

    def _calc_prob(self, z0: np.ndarray, z1: np.ndarray, z2: np.ndarray, ddz1: np.ndarray, dt: float) -> np.ndarray:
        t0 = (b := (z2 - z0) / (2 * dt)) / (2 * (a := ddz1 / 2))
        veff = np.sqrt(np.maximum((g_min := a * t0**2 - b * t0 + z1) * ddz1, 0))
        return np.exp(-0.5 * np.pi * (g_min**2) / np.maximum(veff, 1e-14))

    def _normalize_probs(self, probs: np.ndarray) -> np.ndarray:
        if np.any(mask_gt1 := (row_sums := np.sum(probs, axis=1)) > 1):
            probs[mask_gt1] /= row_sums[mask_gt1][:, np.newaxis]
        return probs

    def _update_history(self, ensemble: Ensemble, H: Hamiltonian, time: float) -> bool:
        self._V_history.append(H.pot.eval_a(list(ensemble.r.T), time))
        return len(self._V_history) == 3

    def jump(self, ensemble: Ensemble, H: Hamiltonian, dt: float, time: float) -> None:
        if not self._update_history(ensemble, H, time): return

        probs = np.zeros((ensemble.ntraj, H.pot.nstate))

        for j in range(H.pot.nstate):
            z0, z1, z2, dz0, dz1, ddz1 = self._calc_gaps(ensemble.states, j, dt)

            if np.any(mask := (dz0 * dz1 <= 0) & (ddz1 > 0) & (ensemble.states != j)):
                probs[mask, j] = self._calc_prob(z0[mask], z1[mask], z2[mask], ddz1[mask], dt)

        probs = self._normalize_probs(probs)
        ensemble.states, ensemble.p = self._apply_jump(ensemble, H, time, probs)
