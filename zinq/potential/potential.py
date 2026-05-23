"""Module defining the abstract base class for potential energy operators."""

from abc import ABC, abstractmethod
from typing import Annotated

import numpy as np
from pydantic import Field


class Potential(ABC):
    """
    Abstract base class for potential energy operators.

    Defines the spatial dimension, number of electronic states, and time dependence,
    as well as diabatic/adiabatic evaluations.
    """

    @property
    @abstractmethod
    def ndim(self) -> int:
        """
        Get the spatial dimensionality of the potential.

        Returns
        -------
        int
            Number of spatial dimensions.

        """

    @property
    @abstractmethod
    def nstate(self) -> int:
        """
        Get the number of electronic states.

        Returns
        -------
        int
            Number of states.

        """

    @property
    def is_td(self) -> bool:
        """
        Check if the potential is time-dependent.

        Returns
        -------
        bool
            True if time-dependent, False otherwise.

        """
        return False

    @abstractmethod
    def eval_d(self, r: list[np.ndarray], time: float) -> np.ndarray:
        """
        Evaluate the potential matrix in the diabatic representation.

        Parameters
        ----------
        r : list[np.ndarray]
            List of coordinate arrays.
        time : float
            Current simulation time.

        Returns
        -------
        np.ndarray
            The diabatic potential energy matrix.

        """

    def eval_a(self, r: list[np.ndarray], time: float) -> np.ndarray:
        """
        Evaluate the potential eigenvalues in the adiabatic representation.

        Parameters
        ----------
        r : list[np.ndarray]
            List of coordinate arrays.
        time : float
            Current simulation time.

        Returns
        -------
        np.ndarray
            The adiabatic potential energy eigenvalues.

        """
        return np.linalg.eigvalsh(self.eval_d(r, time))


from .harmonic import Harmonic  # noqa: E402
from .time_linear import TimeLinear  # noqa: E402
from .tully import TullyFirst  # noqa: E402

AnyPotential = Annotated[Harmonic | TimeLinear | TullyFirst, Field(discriminator="name")]
