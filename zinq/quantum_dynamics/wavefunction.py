"""Module representing the quantum mechanical wavefunction/wavepacket."""

from functools import cached_property

import numpy as np

from .grid import Grid
from .hamiltonian import Hamiltonian
from .initial_conditions import InitialConditions
from .options import InitialConditionsConfig

NORM_THRESHOLD = 1e-14


class Wavefunction:
    """
    Representation of the system wavefunction on a spatial grid.

    Attributes
    ----------
    data : np.ndarray
        The complex wavefunction values on the grid.

    """

    data: np.ndarray

    def __init__(self, *, ic: InitialConditions, grid: Grid, ham: Hamiltonian) -> None:
        """
        Initialize a Gaussian wavepacket according to the initial conditions.

        Parameters
        ----------
        ic : InitialConditions
            The initial wavepacket parameters.
        grid : Grid
            The spatial grid.
        ham : Hamiltonian
            The system Hamiltonian.

        """
        self.data = np.zeros((*[grid.npoint] * grid.ndim, ham.pot.nstate), dtype=np.complex128)

        exponent = np.zeros_like(grid.pos[0], dtype=np.complex128)

        for pos, r, k, g in zip(grid.pos, ic.pos, ic.mom, ic.gamma, strict=False):
            exponent += -0.5 * g * (dr := pos - r) ** 2 + 1j * k * dr

        self.data[..., ic.state] = np.exp(exponent)

        self.normalize(grid)

        if ic.adia:
            self.data = self.to_dia(ham).data

    @classmethod
    def from_data(cls, data: np.ndarray) -> "Wavefunction":
        """
        Create a Wavefunction instance directly from a data array.

        Parameters
        ----------
        data : np.ndarray
            The wavefunction data array.

        Returns
        -------
        Wavefunction
            The initialized Wavefunction instance.

        """
        return (wfn := cls.__new__(cls)).__dict__.update({"data": data}) or wfn

    @classmethod
    def from_options(
        cls, opt: InitialConditionsConfig, grid: Grid, ham: Hamiltonian
    ) -> "Wavefunction":
        """
        Create a Wavefunction instance from InitialConditionsConfig options.

        Parameters
        ----------
        opt : InitialConditionsConfig
            The initial conditions configuration.
        grid : Grid
            The spatial grid.
        ham : Hamiltonian
            The system Hamiltonian.

        Returns
        -------
        Wavefunction
            The initialized Wavefunction instance.

        """
        return cls(ic=InitialConditions.from_options(opt), grid=grid, ham=ham)

    @cached_property
    def ndim(self) -> int:
        """
        Get the spatial dimensionality of the wavefunction.

        Returns
        -------
        int
            Number of spatial dimensions.

        """
        return len(self.data.shape) - 1

    @cached_property
    def nstate(self) -> int:
        """
        Get the number of electronic states.

        Returns
        -------
        int
            Number of electronic states.

        """
        return self.data.shape[-1]

    @property
    def density(self) -> np.ndarray:
        """
        Get the probability density of the wavefunction.

        Returns
        -------
        np.ndarray
            The absolute square of the wavefunction.

        """
        return np.abs(self.data) ** 2

    @property
    def rho(self) -> np.ndarray:
        """
        Get the total spatial probability density (summed over states).

        Returns
        -------
        np.ndarray
            Total spatial probability density.

        """
        return np.sum(self.density, axis=-1)

    def copy(self) -> "Wavefunction":
        """
        Create a deep copy of the wavefunction.

        Returns
        -------
        Wavefunction
            A copy of the wavefunction with copied data.

        """
        return Wavefunction.from_data(self.data.copy())

    def ke(self, grid: Grid, ham: Hamiltonian) -> float:
        """
        Calculate the kinetic energy expectation value.

        Parameters
        ----------
        grid : Grid
            The spatial grid.
        ham : Hamiltonian
            The system Hamiltonian.

        Returns
        -------
        float
            Kinetic energy expectation value.

        """
        return np.sum(self.to_kspace().rho * ham.T) * grid.measure

    def mom(self, grid: Grid) -> np.ndarray:
        """
        Calculate the momentum expectation value vector.

        Parameters
        ----------
        grid : Grid
            The spatial grid.

        Returns
        -------
        np.ndarray
            Momentum expectation values along each dimension.

        """
        rho_k = self.to_kspace().rho

        return np.array([np.sum(rho_k * k) for k in grid.mom]) * grid.measure

    def norm(self, grid: Grid) -> float:
        """
        Calculate the total L2 norm of the wavefunction.

        Parameters
        ----------
        grid : Grid
            The spatial grid.

        Returns
        -------
        float
            The L2 norm of the wavefunction.

        """
        return np.real(np.vdot(self.data, self.data)).item() * grid.measure

    def normalize(self, grid: Grid) -> None:
        """
        Normalize the wavefunction to unity norm.

        Parameters
        ----------
        grid : Grid
            The spatial grid.

        """
        self.data /= np.sqrt(norm) if (norm := self.norm(grid)) > NORM_THRESHOLD else 1

    def overlap(self, other: "Wavefunction", grid: Grid) -> complex:
        """
        Calculate the overlap integral with another wavefunction.

        Parameters
        ----------
        other : Wavefunction
            The other wavefunction.
        grid : Grid
            The spatial grid.

        Returns
        -------
        complex
            The complex overlap integral.

        """
        return np.vdot(self.data, other.data).item() * grid.measure

    def pe(self, grid: Grid, ham: Hamiltonian) -> float:
        """
        Calculate the potential energy expectation value.

        Parameters
        ----------
        grid : Grid
            The spatial grid.
        ham : Hamiltonian
            The system Hamiltonian.

        Returns
        -------
        float
            Potential energy expectation value.

        """
        pe_dens = np.einsum("...i,...ij,...j->...", np.conj(self.data), ham.V, self.data)

        return np.real(np.sum(pe_dens)) * grid.measure

    def pop(self, grid: Grid) -> np.ndarray:
        """
        Calculate the state population probabilities.

        Parameters
        ----------
        grid : Grid
            The spatial grid.

        Returns
        -------
        np.ndarray
            State population probabilities.

        """
        return np.sum(self.density, axis=tuple(range(grid.ndim))) * grid.measure

    def pos(self, grid: Grid) -> np.ndarray:
        """
        Calculate the coordinate expectation value vector.

        Parameters
        ----------
        grid : Grid
            The spatial grid.

        Returns
        -------
        np.ndarray
            Coordinate expectation values along each dimension.

        """
        rho = self.rho

        return np.array([np.sum(rho * r) for r in grid.pos]) * grid.measure

    def project_out(self, others: list["Wavefunction"], grid: Grid) -> None:
        """
        Project out a list of other wavefunctions from this wavefunction.

        Parameters
        ----------
        others : list[Wavefunction]
            The wavefunctions to project out.
        grid : Grid
            The spatial grid.

        """
        for other in others:
            self.data -= np.conj(self.overlap(other, grid)) * other.data

        self.normalize(grid)

    def to_adia(self, ham: Hamiltonian) -> "Wavefunction":
        """
        Transform the wavefunction to the adiabatic basis.

        Parameters
        ----------
        ham : Hamiltonian
            The system Hamiltonian containing the transformation matrix U.

        Returns
        -------
        Wavefunction
            The wavefunction in the adiabatic basis.

        """
        return Wavefunction.from_data(np.einsum("...ji,...j->...i", np.conj(ham.U), self.data))

    def to_dia(self, ham: Hamiltonian) -> "Wavefunction":
        """
        Transform the wavefunction to the diabatic basis.

        Parameters
        ----------
        ham : Hamiltonian
            The system Hamiltonian containing the transformation matrix U.

        Returns
        -------
        Wavefunction
            The wavefunction in the diabatic basis.

        """
        return Wavefunction.from_data(np.einsum("...ij,...j->...i", ham.U, self.data))

    def to_kspace(self) -> "Wavefunction":
        """
        Transform the wavefunction to momentum space (k-space).

        Returns
        -------
        Wavefunction
            The wavefunction in momentum space.

        """
        return Wavefunction.from_data(np.fft.fftn(self.data, axes=range(self.ndim), norm="ortho"))
