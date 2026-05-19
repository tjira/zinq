import numpy as np

from .grid import Grid
from .hamiltonian import Hamiltonian
from .wavefunction import Wavefunction
from ..potential import Potential


class StrangSplit:
    R: np.ndarray
    K: np.ndarray
    unit: complex
    dt: float

    def __init__(self, grid: Grid, H: Hamiltonian, pot: Potential, dt: float, imaginary: bool):
        self.dt, self.unit = dt, -0.5 * (1 if imaginary else 1j) * dt
        self.K = np.exp(2 * self.unit * H.T)[..., None]

        self._update_R(grid, H, pot, 0)

    def step(self, wfn: Wavefunction, grid: Grid, H: Hamiltonian, pot: Potential, time: float):
        if pot.is_time_dependent:
            self._update_R(grid, H, pot, time)

        self._apply_R(wfn)

        wfn.data = np.fft.fftn(wfn.data, axes=range(grid.ndim))
        wfn.data *= self.K
        wfn.data = np.fft.ifftn(wfn.data, axes=range(grid.ndim))

        self._apply_R(wfn)

        if self.unit.imag == 0: wfn.normalize(grid)
    
    def _apply_R(self, wfn):
        wfn.data = np.einsum("...ij,...j->...i", self.R, wfn.data)

    def _update_R(self, grid: Grid, H: Hamiltonian, pot: Potential, time: float):
        H._update_V(grid, pot, time)

        self.R = H.U @ (np.exp(self.unit * H.W)[..., None] * H.U.conj().mT)
