from ...backend import np
from ..ensemble import Ensemble
from ..hamiltonian import Hamiltonian
from .surface_hopping import SurfaceHopping


class LandauZener(SurfaceHopping):
    _V_history: list[np.ndarray]
    _rng: np.random.Generator

    def __init__(self, seed: int) -> None:
        self._V_history = []
        self._rng = np.random.default_rng(seed)

    def jump(self, ensemble: Ensemble, H: Hamiltonian, dt: float, time: float) -> None:
        V = H.pot.eval_a([ensemble.r[:, i] for i in range(ensemble.ndim)], time)

        self._V_history.append(V)

        if len(self._V_history) > 3: self._V_history.pop(0)

        if len(self._V_history) < 3: return

        V0, V1, V2 = self._V_history[2], self._V_history[1], self._V_history[0]

        states, indices, probs = ensemble.states, np.arange(ensemble.ntraj), np.zeros((ensemble.ntraj, H.pot.nstate))

        V0_s, V1_s, V2_s = V0[states, indices], V1[states, indices], V2[states, indices]

        for j in range(H.pot.nstate):

            z0, z1, z2 = np.abs(V0_s - V0[j, :]), np.abs(V1_s - V1[j, :]), np.abs(V2_s - V2[j, :])

            dz0, dz1, ddz1 = (z0 - z1) / dt, (z1 - z2) / dt, (z0 - 2 * z1 + z2) / (dt * dt)

            mask = (dz0 * dz1 <= 0) & (ddz1 > 0) & (states != j)

            if np.any(mask):

                a_coef, b_coef = ddz1[mask] / 2, (z0[mask] - z2[mask]) / (2 * dt)
                t0 = -b_coef / (2 * a_coef)
                g_min = a_coef * t0**2 + b_coef * t0 + z1[mask]

                veff = np.sqrt(np.maximum(g_min * ddz1[mask], 0))
                delta = 0.25 * (g_min**2) / np.maximum(veff, 1e-14)

                probs[mask, j] = np.exp(-2 * np.pi * delta)

        row_sums = np.sum(probs, axis=1)
        mask_gt1 = row_sums > 1
        probs[mask_gt1] /= row_sums[mask_gt1][:, np.newaxis]

        random_vals, cum_probs = self._rng.random(ensemble.ntraj), np.cumsum(probs, axis=1)

        jumps = random_vals[:, np.newaxis] < cum_probs
        jump_happened = np.any(jumps, axis=1)

        if np.any(jump_happened):
            trajs = np.where(jump_happened)[0]
            new_states = np.argmax(jumps, axis=1)[jump_happened]
            
            for i, traj_idx in enumerate(trajs):
                old_state, new_state = ensemble.states[traj_idx], new_states[i]
                
                dV = V0[new_state, traj_idx] - V0[old_state, traj_idx]
                p_sq = np.sum(ensemble.p[traj_idx]**2)
                
                new_p_sq = p_sq - 2 * H.mass * dV
                
                if new_p_sq > 0:
                    ensemble.p[traj_idx] *= np.sqrt(new_p_sq / p_sq)
                    ensemble.states[traj_idx] = new_state
