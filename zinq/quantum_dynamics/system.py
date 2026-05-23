"""Module representing the physical system undergoing quantum dynamics simulation."""

from dataclasses import dataclass

import numpy as np

from .grid import Grid
from .hamiltonian import Hamiltonian
from .options import Options
from .wavefunction import Wavefunction


@dataclass
class System:
    """
    Representation of the entire physical system (wavepacket, Hamiltonian, and grid).

    Attributes
    ----------
    ham : Hamiltonian
        The Hamiltonian operator for the system.
    grid : Grid
        The spatial grid on which the system is defined.
    wfn : Wavefunction
        The wavefunction/wavepacket of the system.

    """

    ham: Hamiltonian
    grid: Grid
    wfn: Wavefunction

    @classmethod
    def from_options(cls, opt: Options) -> "System":
        """
        Create a System instance from simulation options.

        Parameters
        ----------
        opt : Options
            The simulation options.

        Returns
        -------
        System
            The initialized System instance.

        """
        grid = Grid.from_options(opt.grid)
        ham = Hamiltonian.from_options(opt.hamiltonian, grid)
        wfn = Wavefunction.from_options(opt.initial_conditions, grid, ham)
        return cls(grid=grid, ham=ham, wfn=wfn)

    @property
    def ke(self) -> float:
        """
        Get the kinetic energy expectation value.

        Returns
        -------
        float
            Kinetic energy expectation value.

        """
        return self.wfn.ke(self.grid, self.ham)

    @property
    def mom(self) -> np.ndarray:
        """
        Get the momentum expectation value vector.

        Returns
        -------
        np.ndarray
            Momentum expectation values along each dimension.

        """
        return self.wfn.mom(self.grid)

    @property
    def norm(self) -> float:
        """
        Get the norm of the wavefunction.

        Returns
        -------
        float
            The L2 norm of the wavepacket.

        """
        return self.wfn.norm(self.grid)

    @property
    def pe(self) -> float:
        """
        Get the potential energy expectation value.

        Returns
        -------
        float
            Potential energy expectation value.

        """
        return self.wfn.pe(self.grid, self.ham)

    @property
    def pop(self) -> np.ndarray:
        """
        Get the populations of the electronic states.

        Returns
        -------
        np.ndarray
            State population probabilities.

        """
        return self.wfn.pop(self.grid)

    @property
    def pos(self) -> np.ndarray:
        """
        Get the coordinate expectation value vector.

        Returns
        -------
        np.ndarray
            Coordinate expectation values along each dimension.

        """
        return self.wfn.pos(self.grid)

    def normalize(self) -> None:
        """Normalize the wavefunction of the system."""
        self.wfn.normalize(self.grid)

    def update_v(self, time: float) -> None:
        """
        Update the potential energy of the Hamiltonian to the given time.

        Parameters
        ----------
        time : float
            The current simulation time.

        """
        self.ham.update_v(self.grid, time)
