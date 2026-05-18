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

    def _propagate_coefficients(self, ensemble: Ensemble, H_eff: np.ndarray, dot_d: np.ndarray, dt: float) -> None:
        h = dt / self._substeps
        M = -dot_d.astype(complex)
        M[:, np.arange(H_eff.shape[1]), np.arange(H_eff.shape[1])] -= 1j * H_eff

        hM = h * M
        M2 = hM @ hM
        M3 = hM @ M2 / 3
        M4 = hM @ M3 / 4
        m_step = np.eye(H_eff.shape[1]) + hM + M2 / 2 + M3 + M4

        prop = np.linalg.matrix_power(m_step, self._substeps)
        ensemble.c = (prop @ ensemble.c[..., np.newaxis]).squeeze(-1)

    def _calc_probs(self, ensemble: Ensemble, dot_d: np.ndarray, dt: float) -> np.ndarray:
        trajs = np.arange(ensemble.ntraj)
        c_i = ensemble.c[trajs, ensemble.states, np.newaxis]
        term = ensemble.c * np.conj(c_i) * dot_d[trajs, ensemble.states, :]
        probs = 2 * np.real(term) * dt / (np.abs(c_i) ** 2 + 1e-14)

        probs = np.maximum(probs, 0)
        probs[trajs, ensemble.states] = 0

        return self._normalize_probs(probs)

    def _normalize_probs(self, probs: np.ndarray) -> np.ndarray:
        if np.any(mask_gt1 := (row_sums := np.sum(probs, axis=1)) > 1):
            probs[mask_gt1] /= row_sums[mask_gt1][:, np.newaxis]
        return probs

    def jump(self, ensemble: Ensemble, pot: Potential, dt: float, time: float) -> None:
        V_a, U = np.linalg.eigh(pot.eval_d(list(ensemble.r.T), time))
        self._U_history.append(U)

        if len(self._U_history) < 2: return

        dot_d = self._tdc.evaluate(self._U_history[0], self._U_history[1], dt)
        H_eff = V_a - np.mean(V_a, axis=1, keepdims=True)

        self._propagate_coefficients(ensemble, H_eff, dot_d, dt)

        probs = self._calc_probs(ensemble, dot_d, dt)
        ensemble.states, ensemble.p = self._apply_jump(ensemble, V_a, probs)
