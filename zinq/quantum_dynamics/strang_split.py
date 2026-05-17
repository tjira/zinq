from typing import Optional

from ..backend import np
from .grid import Grid
from .hamiltonian import Hamiltonian


class StrangSplit:
    H: Hamiltonian
    R: np.ndarray
    K: np.ndarray
    unit: complex

    def __init__(self, grid: Grid, H: Hamiltonian, dt: float, imaginary: bool):
        assert dt > 0, f"TIME STEP MUST BE POSITIVE, GOT {dt}"

        self.unit, self.H = -0.5 * (1.0 if imaginary else 1j) * dt, H

        self._recalculate_R(grid, 0)
        k_sq = sum((k**2 for k in grid.momentum), np.zeros_like(grid.momentum[0]))
        self.K = np.exp(self.unit * k_sq / H.mass)[..., np.newaxis]

    def _recalculate_R(self, grid: Grid, time: float):
        V = np.moveaxis(self.H.pot.eval_d(grid.position, time), [0, 1], [-2, -1])
        W, U = np.linalg.eigh(V)

        if self.H.cap: W = W - 1j * self.H.eval_cap(grid.position)[..., np.newaxis]

        self.R = U @ (np.exp(self.unit * W)[..., None] * U.conj().mT)

    def step(self, wfn, grid: Grid, time: float) -> np.ndarray:
        if self.H.pot.is_time_dependent:
            self._recalculate_R(grid, time + abs(self.unit.real + self.unit.imag))

        decay = self._apply_R(wfn)
        wfn.data = np.fft.fftn(wfn.data, axes=range(wfn.ndim))
        wfn.data *= self.K
        wfn.data = np.fft.ifftn(wfn.data, axes=range(wfn.ndim))
        decay += self._apply_R(wfn)

        if self.unit.imag == 0: wfn.normalize()

        return decay

    def _apply_R(self, wfn) -> np.ndarray:
        wfn.data = np.einsum("...ij,...j->...i", self.R, wfn.data)

        if not self.H.cap or not self.H.cap.track_population: return np.zeros(wfn.nstate)

        col_norm = np.sum(np.abs(self.R[..., :, 0])**2, axis=-1)

        if np.allclose(col_norm, 1): return np.zeros(wfn.nstate)

        decay_density = (1 / np.maximum(col_norm[..., np.newaxis], 1e-14) - 1) * np.abs(wfn.data)**2

        return np.sum(decay_density, axis=tuple(range(wfn.ndim))) * wfn.measure
