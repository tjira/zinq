"""Module containing the Result class for quantum dynamics simulations."""

from dataclasses import dataclass

import numpy as np


@dataclass(frozen=True, kw_only=True)
class Result:
    """
    Result of a quantum dynamics simulation.

    Attributes
    ----------
    total_energy : np.ndarray
        Total energy at each time step.
    kinetic_energy : np.ndarray
        Kinetic energy at each time step.
    potential_energy : np.ndarray
        Potential energy at each time step.
    position : np.ndarray
        Coordinate expectation value vector at each time step.
    momentum : np.ndarray
        Momentum expectation value vector at each time step.
    norm : np.ndarray
        Wavepacket norm at each time step.
    population : np.ndarray
        Electronic state populations at each time step.

    """

    total_energy: np.ndarray
    kinetic_energy: np.ndarray
    potential_energy: np.ndarray
    position: np.ndarray
    momentum: np.ndarray
    norm: np.ndarray
    population: np.ndarray
