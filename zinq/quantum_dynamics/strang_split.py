from ..backend import np
from .grid import Grid
from .hamiltonian import Hamiltonian


class StrangSplit:
    H: Hamiltonian
    R: np.ndarray
    K: np.ndarray
    unit: complex

    U: np.ndarray = None

    def __init__(self, grid: Grid, ham: Hamiltonian, dt: float, imaginary: bool, adiabatic: bool):
        assert dt > 0, f"TIME STEP MUST BE POSITIVE, GOT {dt}"

        self.unit, self.H = -0.5 * (1.0 if imaginary else 1j) * dt, ham

        V = np.moveaxis(ham.potential.eval_d(grid.position), [0, 1], [-2, -1])

        k_squared, W, U = sum(k**2 for k in grid.momentum), *np.linalg.eigh(V)

        if adiabatic: self.U = U

        if ham.cap: W = W - 1j * ham.eval_cap(grid.position)[..., np.newaxis]

        self.R = U @ (np.exp(self.unit * W)[..., None] * U.conj().mT)
        self.K = np.exp(self.unit * k_squared / ham.mass)[..., np.newaxis]

    def step(self, wfn) -> np.ndarray:
        decay = self._apply_R(wfn)

        wfn.data = np.fft.fftn(wfn.data, axes=range(wfn.ndim))
        wfn.data *= self.K
        wfn.data = np.fft.ifftn(wfn.data, axes=range(wfn.ndim))

        decay += self._apply_R(wfn)

        if self.unit.imag == 0: wfn.normalize()

        return decay

    def _apply_R(self, wfn) -> np.ndarray:
        wfn.data = np.einsum("...ij,...j->...i", self.R, wfn.data)

        if not self.H.cap or not self.H.cap.track_population: return 0

        col_norm = np.sum(np.abs(self.R[..., :, 0])**2, axis=-1)

        if np.allclose(col_norm, 1): return 0

        data = np.einsum("...ji,...j->...i", self.U.conj(), wfn.data) if self.U is not None else wfn.data

        return np.sum((1 / col_norm[..., np.newaxis] - 1) * np.abs(data)**2, axis=tuple(range(wfn.ndim))) * wfn.measure

