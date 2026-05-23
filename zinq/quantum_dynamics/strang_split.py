from functools import cached_property

import numpy as np

from .hamiltonian import Hamiltonian
from .system import System
from .wavefunction import Wavefunction


class StrangSplit:
    def __init__(self, *, ham: Hamiltonian, dt: float, imag: bool, adia: bool = False):
        self.dt, self.unit, self.adia = dt, -0.5 * (1 if imag else 1j) * dt, adia

        self._update_R(ham)
        self._update_K(ham)

    @cached_property
    def imag(self):
        return self.unit.imag == 0

    def step(self, system: System, time: float) -> np.ndarray:
        if system.ham.pot.is_td: system.update_V(time); self._update_R(system.ham)

        decay = self._apply_R(system) + (self._apply_K(system.wfn) or 0) + self._apply_R(system)

        if self.imag: system.normalize()

        return decay

    def _apply_K(self, wfn: Wavefunction):
        wfn.data = np.fft.ifftn(np.fft.fftn(wfn.data, axes=range(wfn.ndim)) * self.K, axes=range(wfn.ndim))

    def _apply_R(self, system: System) -> np.ndarray:
        system.wfn.data = np.einsum("...ij,...j->...i", self.R, system.wfn.data)

        if self.imag or not system.ham.absorber: return np.zeros(system.wfn.nstate)

        return self._get_decay(system)

    def _get_decay(self, system: System) -> np.ndarray:
        wfn = system.wfn.to_adia(system.ham) if self.adia else system.wfn

        remaining_norm = np.sum(np.abs(self.R[..., :, 0])**2, axis=-1)

        if np.allclose(remaining_norm, 1): return np.zeros(wfn.nstate)

        decay_density = (1 / np.maximum(remaining_norm[..., np.newaxis], 1e-14) - 1) * wfn.density

        return np.sum(decay_density, axis=tuple(range(wfn.ndim))) * system.grid.measure

    def _update_K(self, ham: Hamiltonian):
        self.K = np.exp(2 * self.unit * ham.T)[..., None]

    def _update_R(self, ham: Hamiltonian):
        self.R = ham.U @ (np.exp(self.unit * ham.W)[..., np.newaxis] * ham.U.conj().mT)
