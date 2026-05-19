import numpy as np

from .grid import Grid
from .initial_conditions import InitialConditions


class Wavefunction:
    data: np.ndarray

    def __init__(self, ic: InitialConditions, grid: Grid, nstate: int):
        self.data = np.empty((grid.npoint**grid.ndim, nstate), dtype=np.complex128)

        self.data.fill(0)

        exponent = np.zeros_like(grid.position[0], dtype=np.complex128)

        for i in range(grid.ndim):
            dr = grid.position[i] - ic.position[i]
            exponent += -0.5 * ic.gamma[i] * dr**2 + 1j * ic.momentum[i] * dr

        self.data[..., ic.state] = np.exp(exponent)

        # self.normalize()
