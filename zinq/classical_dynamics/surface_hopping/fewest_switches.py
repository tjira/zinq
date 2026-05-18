from collections import deque

from ...backend import np
from ..ensemble import Ensemble
from ..time_derivative_coupling import TimeDerivativeCoupling
from .surface_hopping import SurfaceHopping
from ...potential import Potential


class FewestSwitches(SurfaceHopping):
    _tdc: TimeDerivativeCoupling
    _substeps: int
    _U_history: deque

    def __init__(self, seed: int, tdc: TimeDerivativeCoupling, substeps: int = 10) -> None:
        super().__init__(seed)
        self._tdc = tdc
        self._substeps = substeps
        self._U_history = deque(maxlen=2)

    def _coef_propagator(self, H_eff: np.ndarray, dot_d: np.ndarray, dt: float) -> np.ndarray:
        M = -dot_d - 1j * H_eff[..., np.newaxis] * np.eye(H_eff.shape[1])

        M1 = dt / self._substeps * M
        M2 = M1 @ M1
        M3 = M1 @ M2
        M4 = M1 @ M3

        return np.eye(H_eff.shape[1]) + M1 + M2 / 2 + M3 / 6 + M4 / 24

    def _calc_probs(self, ensemble: Ensemble, dot_d: np.ndarray, dt: float) -> np.ndarray:
        indices = np.arange(ensemble.ntraj)

        c_active = ensemble.c[indices, ensemble.states, None] 
        coupling = dot_d[indices, ensemble.states, :]

        prob_flux = 2 * np.real(ensemble.c * np.conj(c_active) * coupling) * dt
        probs = prob_flux / (np.abs(c_active)**2 + 1e-14)

        probs = np.maximum(probs, 0)
        probs[indices, ensemble.states] = 0

        return self._normalize_probs(probs)

    def _integrate_substeps(self, ensemble: Ensemble, dot_d: np.ndarray, prop: np.ndarray, dt: float) -> tuple[np.ndarray, np.ndarray]:
        jump_mask, target_states = np.zeros(ensemble.ntraj, dtype=bool), ensemble.states.copy()

        for _ in range(self._substeps):
            ensemble.c = np.einsum("nij,nj->ni", prop, ensemble.c)

            probs = self._calc_probs(ensemble, dot_d, dt / self._substeps)
            step_jump_mask, step_target_states = self._propose_jumps(probs)

            if np.any(valid := step_jump_mask & ~jump_mask):
                target_states[valid], jump_mask = step_target_states[valid], jump_mask | valid

        return jump_mask, target_states

    def _normalize_probs(self, probs: np.ndarray) -> np.ndarray:
        if np.any(mask_gt1 := (row_sums := np.sum(probs, axis=1)) > 1):
            probs[mask_gt1] /= row_sums[mask_gt1][:, np.newaxis]
        return probs

    def _prepare_propagator(self, V_a: np.ndarray, dt: float) -> tuple[np.ndarray, np.ndarray]:
        dot_d = self._tdc.evaluate(self._U_history[0], self._U_history[1], dt)
        
        H_eff = V_a - np.mean(V_a, axis=1, keepdims=True)
        prop = self._coef_propagator(H_eff, dot_d, dt)
        
        return dot_d, prop

    def _update_history(self, ensemble: Ensemble, pot: Potential, time: float) -> np.ndarray:
        V_a, U = np.linalg.eigh(pot.eval_d(list(ensemble.r.T), time))

        if len(self._U_history) > 0:
            U = np.where(np.sum(self._U_history[-1] * U, axis=-2, keepdims=True) < 0, -U, U)

        self._U_history.append(U)

        return V_a

    def jump(self, ensemble: Ensemble, pot: Potential, dt: float, time: float) -> None:
        V_a = self._update_history(ensemble, pot, time)

        if len(self._U_history) < 2: return

        dot_d, prop = self._prepare_propagator(V_a, dt)

        jump_mask, target_states = self._integrate_substeps(ensemble, dot_d, prop, dt)

        self._apply_jump(ensemble, V_a, jump_mask, target_states)
