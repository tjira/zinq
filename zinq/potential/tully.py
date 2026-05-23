"""Module representing Tully's first simple avoided crossing potential model."""

from dataclasses import dataclass
from typing import Literal

import numpy as np

from .potential import Potential


@dataclass(frozen=True, kw_only=True)
class TullyFirst(Potential):
    """
    Representation of Tully's first avoided crossing potential model.

    Attributes
    ----------
    name : Literal["tully_1"]
        The name of the potential.
    A : float
        Asymptotic potential energy parameter.
    B : float
        Potential transition steepness parameter.
    C : float
        Diabatic coupling strength parameter.
    D : float
        Diabatic coupling width parameter.

    """

    name: Literal["tully_1"] = "tully_1"

    A: float = 0.01
    B: float = 1.6
    C: float = 0.005
    D: float = 1

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

    def eval_d(self, r: list[np.ndarray], time: float = 0.0) -> np.ndarray:  # noqa: ARG002
        """
        Evaluate the diabatic potential energy matrix.

        Parameters
        ----------
        r : list[np.ndarray]
            List of coordinate arrays.
        time : float
            Current simulation time (unused).

        Returns
        -------
        np.ndarray
            The diabatic potential energy matrix.

        """
        v00 = np.sign(r[0]) * self.A * (1 - np.exp(-self.B * np.abs(r[0])))
        v01 = self.C * np.exp(-self.D * r[0] ** 2)
        v11 = -v00

        return np.stack(
            [np.stack([v00, v01], axis=-1), np.stack([v01, v11], axis=-1)],
            axis=-2,
        )
