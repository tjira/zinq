"""Module representing the Hamiltonian operator for quantum dynamics."""

import numpy as np

from zinq.potential import Potential

from .absorber import Absorber
from .grid import Grid
from .options import HamiltonianConfig


class Hamiltonian:
    """
    Representation of the Hamiltonian operator.

    Attributes
    ----------
    pot : Potential
        The potential energy operator.
    m : float
        The particle mass.
    absorber : Absorber | None
        The boundary absorber if configured.
    T : np.ndarray
        The kinetic energy operator in momentum space.
    W : np.ndarray
        Eigenvalues of the potential energy matrix (adiabatic potentials).
    U : np.ndarray
        Eigenvectors of the potential energy matrix (diabatic-to-adiabatic transformation).
    V : np.ndarray
        Potential energy matrix in the diabatic basis.

    """

    pot: Potential
    m: float
    absorber: Absorber | None
    T: np.ndarray
    W: np.ndarray
    U: np.ndarray
    V: np.ndarray

    def __init__(
        self,
        *,
        grid: Grid,
        pot: Potential,
        m: float,
        absorber: Absorber | None = None,
    ) -> None:
        """
        Initialize the Hamiltonian operator.

        Parameters
        ----------
        grid : Grid
            The spatial simulation grid.
        pot : Potential
            The potential energy operator.
        m : float
            The particle mass.
        absorber : Absorber | None
            Optional absorbing boundary potential.

        """
        self.pot, self.m, self.absorber = pot, m, absorber

        self.update_v(grid, 0)

        self.T = sum((k**2 for k in grid.mom), np.zeros_like(grid.mom[0])) / (2 * m)

    @classmethod
    def from_options(cls, opt: HamiltonianConfig, grid: Grid) -> "Hamiltonian":
        """
        Create a Hamiltonian instance from HamiltonianConfig options.

        Parameters
        ----------
        opt : HamiltonianConfig
            The Hamiltonian configuration.
        grid : Grid
            The simulation grid.

        Returns
        -------
        Hamiltonian
            The initialized Hamiltonian operator.

        """
        absorber = Absorber.from_options(opt.absorber) if opt.absorber else None

        return cls(grid=grid, pot=opt.potential, m=opt.mass, absorber=absorber)

    def update_v(self, grid: Grid, time: float) -> None:
        """
        Update the potential energy operator V and its eigensystem.

        Parameters
        ----------
        grid : Grid
            The simulation grid.
        time : float
            The current simulation time.

        """
        v_eval = self.pot.eval_d(grid.pos, time)
        self.W, self.U, self.V = *np.linalg.eigh(v_eval), v_eval

        self._fix_gauge()

        if self.absorber:
            self.W = self.W - 1j * self.absorber.eval(grid.pos)[..., np.newaxis]

    def _fix_gauge(self) -> None:
        """Fix the gauge/phases of the electronic states to ensure smooth continuity."""
        if self.pot.ndim != 1:
            return

        flips = np.where(self._get_overlaps() < 0, -1, 1)
        self.U *= np.cumprod(np.insert(flips, 0, 1, axis=0), axis=0)[:, np.newaxis, :]

    def _get_overlaps(self) -> np.ndarray:
        """Calculate state overlaps between adjacent grid points."""
        return np.einsum("...ij,...ij->...j", self.U[:-1], self.U[1:])
