from collections import deque

from ...backend import np
from ..ensemble import Ensemble
from ..hamiltonian import Hamiltonian
from ..time_derivative_coupling import TimeDerivativeCoupling
from .surface_hopping import SurfaceHopping


class FewestSwitches(SurfaceHopping):
    _tdc: TimeDerivativeCoupling
    _substeps: int
    _U_history: deque

    def __init__(self, seed: int, tdc: TimeDerivativeCoupling, substeps: int = 10) -> None:
        super().__init__(seed)
        self._tdc = tdc
        self._substeps = substeps
        self._U_history = deque(maxlen=2)

    def _deriv(self, c: np.ndarray, V: np.ndarray, V_avg: np.ndarray, dot_d: np.ndarray) -> np.ndarray:
        return -1j * (V - V_avg) * c - np.einsum("nij,nj->ni", dot_d, c)

    def jump(self, ensemble: Ensemble, H: Hamiltonian, dt: float, time: float) -> None:
        V_d = H.pot.eval_d(list(ensemble.r.T), time)
        V_a, U = np.linalg.eigh(V_d)

        self._U_history.append(U)
        if len(self._U_history) < 2:
            return

        dot_d = self._tdc.evaluate(self._U_history[0], self._U_history[1], dt)

        V_avg = np.mean(V_a, axis=1, keepdims=True)
        dt_q = dt / self._substeps

        for _ in range(self._substeps):
            k1 = self._deriv(ensemble.c, V_a, V_avg, dot_d)
            k2 = self._deriv(ensemble.c + 0.5 * dt_q * k1, V_a, V_avg, dot_d)
            k3 = self._deriv(ensemble.c + 0.5 * dt_q * k2, V_a, V_avg, dot_d)
            k4 = self._deriv(ensemble.c + dt_q * k3, V_a, V_avg, dot_d)
            ensemble.c += (dt_q / 6.0) * (k1 + 2 * k2 + 2 * k3 + k4)

        probs = np.zeros((ensemble.ntraj, H.pot.nstate))
        states = ensemble.states
        trajs = np.arange(ensemble.ntraj)

        c_curr = ensemble.c[trajs, states]
        denominator = np.abs(c_curr) ** 2 + 1e-14

        for j in range(H.pot.nstate):
            mask = j != states
            re = np.real(ensemble.c[trajs, j] * np.conj(c_curr))
            p = 2 * dot_d[trajs, states, j] * re / denominator * dt
            probs[:, j] = np.where(mask, np.maximum(p, 0), 0)

        if np.any(mask_gt1 := (row_sums := np.sum(probs, axis=1)) > 1):
            probs[mask_gt1] /= row_sums[mask_gt1][:, np.newaxis]

        ensemble.states, ensemble.p = self._apply_jump(ensemble, H, time, probs)
