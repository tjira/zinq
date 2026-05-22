import numpy as np

from ..potential import Potential
from .absorber import Absorber
from .grid import Grid


class Hamiltonian:
    def __init__(self, grid: Grid, pot: Potential, m: float, absorber: Absorber | None = None):
        self.pot, self.m, self.absorber = pot, m, absorber
        
        self.update_V(grid, 0)

        self.T = sum((k**2 for k in grid.mom), np.zeros_like(grid.mom[0])) / (2 * m)

    def update_V(self, grid: Grid, time: float):
        self.W, self.U, self.V = *np.linalg.eigh(V := self.pot.eval_d(grid.pos, time)), V

        self._fix_gauge()

        if self.absorber:
            self.W = self.W - 1j * self.absorber.eval(grid.pos)[..., np.newaxis]

    def _fix_gauge(self) -> None:
        def get_flips():
            return np.where(self._get_overlaps() < 0, -1, 1)

        if self.pot.ndim != 1: return

        self.U *= np.cumprod(np.insert(get_flips(), 0, 1, axis=0), axis=0)[:, np.newaxis, :]

    def _get_overlaps(self) -> np.ndarray:
        return np.einsum("...ij,...ij->...j", self.U[:-1], self.U[1:])
