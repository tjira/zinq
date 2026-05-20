import numpy as np

from ..potential import Potential
from .grid import Grid
from .hamiltonian import Hamiltonian
from .wavefunction import Wavefunction


class StrangSplit:
    R: np.ndarray
    K: np.ndarray
    unit: complex
    dt: float

    def __init__(self, H: Hamiltonian, dt: float, imaginary: bool):
        self.dt, self.unit = dt, -0.5 * (1 if imaginary else 1j) * dt

        self._update_R(H)
        self._update_K(H)

    def step(self, wfn: Wavefunction, grid: Grid, H: Hamiltonian, pot: Potential):
        if pot.is_time_dependent: self._update_R(H)

        self._apply_R(wfn)
        self._apply_K(wfn, grid.ndim)
        self._apply_R(wfn)

        if self.unit.imag == 0: wfn.normalize(grid)

    def _apply_K(self, wfn, ndim):
        wfn.data = np.fft.ifftn(np.fft.fftn(wfn.data, axes=range(ndim)) * self.K, axes=range(ndim))
    
    def _apply_R(self, wfn):
        wfn.data = np.einsum("...ij,...j->...i", self.R, wfn.data)

    def _update_K(self, H: Hamiltonian):
        self.K = np.exp(2 * self.unit * H.T)[..., None]

    def _update_R(self, H: Hamiltonian):
        self.R = H.U @ (np.exp(self.unit * H.W)[..., None] * H.U.conj().mT)
