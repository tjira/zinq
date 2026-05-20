import numpy as np

from ..potential import Potential
from .grid import Grid


class Hamiltonian:
    pot: Potential
    V: np.ndarray
    W: np.ndarray
    U: np.ndarray
    T: np.ndarray
    m: float

    def __init__(self, grid: Grid, pot: Potential, m: float):
        self.pot, self.m = pot, m
        
        self.update_V(grid, 0)

        self.T = sum((k**2 for k in grid.mom), np.zeros_like(grid.mom[0])) / (2 * m)

    def update_V(self, grid: Grid, time: float):
        self.W, self.U, self.V = *np.linalg.eigh(V := self.pot.eval_d(grid.pos, time)), V
