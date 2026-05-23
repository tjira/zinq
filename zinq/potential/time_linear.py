"""Module representing a time-dependent, linearly varying potential coupling model."""

from dataclasses import dataclass
from typing import Literal

import numpy as np

from .potential import Potential


@dataclass(frozen=True, kw_only=True)
class TimeLinear(Potential):
    """
    Representation of a time-dependent linear coupling potential model.

    Attributes
    ----------
    name : Literal["time_linear"]
        The name of the potential.
    coupling : float
        The constant coupling strength between states.
    slope : float
        The linear rate of potential energy shift over time.

    """

    name: Literal["time_linear"] = "time_linear"

    coupling: float = 2
    slope: float = 10

    @property
    def ndim(self) -> int:
        """
        Get the spatial dimensionality.

        Returns
        -------
        int
            Number of spatial dimensions (always 1).

        """
        return 1

    @property
    def nstate(self) -> int:
        """
        Get the number of electronic states.

        Returns
        -------
        int
            Number of electronic states (always 2).

        """
        return 2

    @property
    def is_td(self) -> bool:
        """
        Check if the potential is time-dependent.

        Returns
        -------
        bool
            True (always time-dependent).

        """
        return True

    def eval_d(self, r: list[np.ndarray], time: float = 0.0) -> np.ndarray:
        """
        Evaluate the diabatic potential energy matrix.

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
        v00 = np.full_like(r[0], self.slope * (time - self.slope))
        v01 = np.full_like(r[0], self.coupling)
        v11 = np.full_like(r[0], -v00)

        return np.stack(
            [np.stack([v00, v01], axis=-1), np.stack([v01, v11], axis=-1)],
            axis=-2,
        )
