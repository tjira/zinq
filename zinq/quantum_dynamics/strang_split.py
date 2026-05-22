from functools import cached_property

import numpy as np

from .grid import Grid
from .hamiltonian import Hamiltonian
from .wavefunction import Wavefunction


class StrangSplit:
    def __init__(self, H: Hamiltonian, dt: float, imag: bool, adiabatic: bool = False):
        self.dt, self.unit, self.adiabatic = dt, -0.5 * (1 if imag else 1j) * dt, adiabatic

        self._update_R(H)
        self._update_K(H)

    @cached_property
    def imag(self):
        return self.unit.imag == 0

    def step(self, wfn: Wavefunction, grid: Grid, H: Hamiltonian, time: float):
        if H.pot.is_td: H.update_V(grid, time); self._update_R(H)

        wfn.absorbed += self._apply_R(wfn, grid, H)
        self._apply_K(wfn)
        wfn.absorbed += self._apply_R(wfn, grid, H)

        if self.imag: wfn.normalize(grid)

    def _apply_K(self, wfn):
        wfn.data = np.fft.ifftn(np.fft.fftn(wfn.data, axes=range(wfn.ndim)) * self.K, axes=range(wfn.ndim))
    
    def _apply_R(self, wfn: Wavefunction, grid: Grid, H: Hamiltonian) -> np.ndarray:
        wfn.data = np.einsum("...ij,...j->...i", self.R, wfn.data)

        if self.imag or not H.absorber: return np.zeros(wfn.nstate)

        return self._get_decay(wfn.to_adia(H) if self.adiabatic else wfn, grid)

    def _get_decay(self, wfn: Wavefunction, grid: Grid) -> np.ndarray:
        remaining_norm = np.sum(np.abs(self.R[..., :, 0])**2, axis=-1)

        if np.allclose(remaining_norm, 1): return np.zeros(wfn.nstate)

        decay_density = (1 / np.maximum(remaining_norm[..., np.newaxis], 1e-14) - 1) * wfn.density

        return np.sum(decay_density, axis=tuple(range(wfn.ndim))) * grid.measure

    def _update_K(self, H: Hamiltonian):
        self.K = np.exp(2 * self.unit * H.T)[..., None]

    def _update_R(self, H: Hamiltonian):
        self.R = H.U @ (np.exp(self.unit * H.W)[..., np.newaxis] * H.U.conj().mT)
