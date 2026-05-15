from ..backend import np
from .grid import Grid
from .hamiltonian import Hamiltonian


class StrangSplit:
    R: np.ndarray
    K: np.ndarray
    unit: complex

    def __init__(self, grid: Grid, ham: Hamiltonian, dt: float, imaginary: bool):
        assert dt > 0, f"TIME STEP MUST BE POSITIVE, GOT {dt}"

        self.unit = -0.5 * (1.0 if imaginary else 1j) * dt

        V = np.moveaxis(ham.potential.eval_d(grid.position), [0, 1], [-2, -1])

        k_squared, W, U = sum(k**2 for k in grid.momentum), *np.linalg.eigh(V)

        if ham.cap:
            W = W - 1j * ham.eval_cap(grid.position)[..., np.newaxis]

        self.R = U @ (np.exp(self.unit * W)[..., None] * U.conj().mT)
        self.K = np.exp(self.unit * k_squared / ham.mass)[..., np.newaxis]

    def step(self, wfn):
        wfn.data = np.einsum("...ij,...j->...i", self.R, wfn.data)
        wfn.data = np.fft.fftn(wfn.data, axes=range(wfn.ndim))
        wfn.data *= self.K
        wfn.data = np.fft.ifftn(wfn.data, axes=range(wfn.ndim))
        wfn.data = np.einsum("...ij,...j->...i", self.R, wfn.data)

        if self.unit.imag == 0: wfn.normalize()
