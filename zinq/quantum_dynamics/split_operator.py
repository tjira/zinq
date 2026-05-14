import numpy as np
import scipy.linalg

from ..potential import Potential


class SplitOperator:
    R: np.ndarray
    K: np.ndarray
    unit: complex

    def __init__(self, position_grid: list[np.ndarray], momentum_grid: list[np.ndarray], potential: Potential, dt: float, mass: float, imaginary: bool = False):
        self.unit = -0.5 * (1.0 if imaginary else 1j) * dt

        self.R = scipy.linalg.expm(self.unit * np.moveaxis(potential.evaluateDiabatic(position_grid), [0, 1], [-2, -1]))
        self.K = np.exp(self.unit * sum(m**2 for m in momentum_grid) / mass)[..., np.newaxis]

    def step(self, wfn):
        wfn.data = np.einsum("...ij,...j->...i", self.R, wfn.data)
        wfn.data = np.fft.fftn(wfn.data, axes=range(wfn.ndim))
        wfn.data *= self.K
        wfn.data = np.fft.ifftn(wfn.data, axes=range(wfn.ndim))
        wfn.data = np.einsum("...ij,...j->...i", self.R, wfn.data)

        if self.unit.imag == 0: wfn.normalize()
