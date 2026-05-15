from dataclasses import dataclass

import numpy as np
import scipy.linalg

from ..potential import Potential


class StrangSplit:
    R: np.ndarray
    K: np.ndarray
    unit: complex

    def __init__(self, /,
        position_grid: list[np.ndarray],
        momentum_grid: list[np.ndarray],
        potential: Potential,
        mass: float,
        dt: float,
        imaginary: bool
    ):
        assert len(position_grid) == potential.ndim, (
            f"POSITION GRID DIMENSION DOES NOT MATCH POTENTIAL DIMENSION"
        )
        assert len(momentum_grid) == potential.ndim, (
            f"MOMENTUM GRID DIMENSION DOES NOT MATCH POTENTIAL DIMENSION"
        )
        assert mass > 0, f"MASS MUST BE POSITIVE, GOT {mass}"
        assert dt > 0, f"TIME STEP MUST BE POSITIVE, GOT {dt}"

        self.unit = -0.5 * (1.0 if imaginary else 1j) * dt

        V = np.moveaxis(potential.evaluate_diabatic(position_grid), [0, 1], [-2, -1])

        k_squared = sum(k**2 for k in momentum_grid)

        self.R = scipy.linalg.expm(self.unit * V)
        self.K = np.exp(self.unit * k_squared / mass)[..., np.newaxis]

    def step(self, wfn):
        wfn.data = np.einsum("...ij,...j->...i", self.R, wfn.data)
        wfn.data = np.fft.fftn(wfn.data, axes=range(wfn.ndim))
        wfn.data *= self.K
        wfn.data = np.fft.ifftn(wfn.data, axes=range(wfn.ndim))
        wfn.data = np.einsum("...ij,...j->...i", self.R, wfn.data)

        if self.unit.imag == 0: wfn.normalize()
