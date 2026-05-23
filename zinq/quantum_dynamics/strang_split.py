"""Module implementing the Strang split-operator propagator for quantum dynamics."""

from functools import cached_property

import numpy as np

from .hamiltonian import Hamiltonian
from .system import System
from .wavefunction import Wavefunction


class StrangSplit:
    """
    Strang split-operator propagator for real/imaginary time wavepacket propagation.

    Attributes
    ----------
    dt : float
        The propagation time step.
    unit : complex
        The scaling factor for propagation (imaginary or real depending on mode).
    adia : bool
        Whether to propagate in the adiabatic representation.
    K : np.ndarray
        The kinetic energy propagator term in momentum space.
    R : np.ndarray
        The potential energy propagator matrix in coordinate space.

    """

    dt: float
    unit: complex
    adia: bool
    K: np.ndarray
    R: np.ndarray

    def __init__(self, *, ham: Hamiltonian, dt: float, imag: bool, adia: bool = False) -> None:
        """
        Initialize the Strang split-operator propagator.

        Parameters
        ----------
        ham : Hamiltonian
            The system Hamiltonian.
        dt : float
            The time step for propagation.
        imag : bool
            True if imaginary time propagation, False for real time propagation.
        adia : bool
            True if propagation should occur in the adiabatic representation.

        """
        self.dt, self.unit, self.adia = dt, -0.5 * (1 if imag else 1j) * dt, adia

        self._update_r(ham)
        self._update_k(ham)

    @cached_property
    def imag(self) -> bool:
        """
        Check if the propagation is in imaginary time.

        Returns
        -------
        bool
            True if propagating in imaginary time, False otherwise.

        """
        return self.unit.imag == 0

    def step(self, system: System, time: float) -> np.ndarray:
        """
        Perform a single Strang split-operator propagation step.

        Parameters
        ----------
        system : System
            The physical system state.
        time : float
            The current simulation time.

        Returns
        -------
        np.ndarray
            The norm decay of electronic states due to the absorber.

        """
        if system.ham.pot.is_td:
            system.update_v(time)
            self._update_r(system.ham)

        decay = self._apply_r(system) + (self._apply_k(system.wfn) or 0) + self._apply_r(system)

        if self.imag:
            system.normalize()

        return decay

    def _apply_k(self, wfn: Wavefunction) -> None:
        """Apply the kinetic energy step in momentum space."""
        wfn.data = np.fft.ifftn(
            np.fft.fftn(wfn.data, axes=range(wfn.ndim)) * self.K, axes=range(wfn.ndim)
        )

    def _apply_r(self, system: System) -> np.ndarray:
        """Apply the potential energy step in coordinate space."""
        system.wfn.data = np.einsum("...ij,...j->...i", self.R, system.wfn.data)

        if self.imag or not system.ham.absorber:
            return np.zeros(system.wfn.nstate)

        return self._get_decay(system)

    def _get_decay(self, system: System) -> np.ndarray:
        """Compute the norm decay per electronic state caused by the absorbing boundary."""
        wfn = system.wfn.to_adia(system.ham) if self.adia else system.wfn

        remaining_norm = np.sum(np.abs(self.R[..., :, 0]) ** 2, axis=-1)

        if np.allclose(remaining_norm, 1):
            return np.zeros(wfn.nstate)

        decay_density = (1 / np.maximum(remaining_norm[..., np.newaxis], 1e-14) - 1) * wfn.density

        return np.sum(decay_density, axis=tuple(range(wfn.ndim))) * system.grid.measure

    def _update_k(self, ham: Hamiltonian) -> None:
        """Update the kinetic energy propagator operator."""
        self.K = np.exp(2 * self.unit * ham.T)[..., None]

    def _update_r(self, ham: Hamiltonian) -> None:
        """Update the potential energy propagator operator matrix."""
        self.R = ham.U @ (np.exp(self.unit * ham.W)[..., np.newaxis] * ham.U.conj().mT)
