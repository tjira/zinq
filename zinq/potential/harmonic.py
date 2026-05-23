"""Module representing a multi-dimensional harmonic oscillator potential."""

from dataclasses import dataclass, field
from typing import Literal

import numpy as np

from .potential import Potential


@dataclass(frozen=True, kw_only=True)
class Harmonic(Potential):
    """
    Representation of a multi-dimensional harmonic potential energy surface.

    Attributes
    ----------
    name : Literal["harmonic"]
        The name of the potential.
    k : list[float]
        Force constants along each dimension.

    """

    name: Literal["harmonic"] = "harmonic"

    k: list[float] = field(default_factory=lambda: [1])

    @property
    def ndim(self) -> int:
        """
        Get the spatial dimensionality.

        Returns
        -------
        int
            Number of spatial dimensions.

        """
        return len(self.k)

    @property
    def nstate(self) -> int:
        """
        Get the number of electronic states.

        Returns
        -------
        int
            Number of electronic states.

        """
        return 1

    def eval_d(self, r: list[np.ndarray], time: float = 0.0) -> np.ndarray:  # noqa: ARG002
        """
        Evaluate the harmonic potential in the diabatic basis.

        Parameters
        ----------
        r : list[np.ndarray]
            List of coordinate arrays.
        time : float
            Current simulation time (unused).

        Returns
        -------
        np.ndarray
            The potential energy matrix.

        """
        v_val = sum(
            (0.5 * self.k[i] * r[i] ** 2 for i in range(self.ndim)),
            start=np.zeros_like(r[0]),
        )

        return v_val[..., np.newaxis, np.newaxis]
