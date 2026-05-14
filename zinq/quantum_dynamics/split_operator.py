import numpy as np
import scipy.linalg

class SplitOperator:
    dt: float
    mass: float
    unit: complex

    def __init__(self, dt: float, mass: float, imaginary: bool = False):
        self.dt, self.mass, self.unit = dt, mass, -0.5 * (1.0 if imaginary else 1j) * dt

    def step(self, wfn, position_grid: list[np.ndarray], momentum_grid: list[np.ndarray], potential):
        V = potential.evaluateDiabatic(position_grid).transpose(*range(2, wfn.data.ndim + 1), 0, 1)
        R = np.array([scipy.linalg.expm(self.unit * v) for v in V.reshape(-1, wfn.nstate, wfn.nstate)]).reshape(*wfn.data.shape[:-1], wfn.nstate, wfn.nstate)
        K = np.exp(self.unit * sum(m**2 for m in momentum_grid) / self.mass)[..., np.newaxis]

        wfn.data = np.einsum("...ij,...j->...i", R, wfn.data)
        wfn.data = np.fft.fftn(wfn.data, axes=range(wfn.ndim))
        wfn.data *= K
        wfn.data = np.fft.ifftn(wfn.data, axes=range(wfn.ndim))
        wfn.data = np.einsum("...ij,...j->...i", R, wfn.data)

        if self.unit.imag == 0: wfn.normalize()
