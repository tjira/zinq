"""Module defining the boundary absorbing potential."""

from dataclasses import dataclass

import numpy as np

from .options import AbsorberConfig


@dataclass(frozen=True, kw_only=True)
class Absorber:
    """
    Representation of an absorbing boundary potential.

    Attributes
    ----------
    limits : np.ndarray
        Spatial boundaries beyond which wavepacket is absorbed. Shape (ndim, 2).
    exponent : float
        Exponential strength factor for absorption.

    """

    limits: np.ndarray
    exponent: float

    @classmethod
    def from_options(cls, opt: AbsorberConfig) -> "Absorber":
        """
        Create an Absorber instance from configuration options.

        Parameters
        ----------
        opt : AbsorberConfig
            The absorber configuration.

        Returns
        -------
        Absorber
            The initialized Absorber instance.

        """
        return cls(
            limits=np.array(opt.limits),
            exponent=opt.exponent,
        )

    def eval(self, r: list[np.ndarray]) -> np.ndarray:
        """
        Evaluate the absorbing potential penalty on the coordinate grids.

        Parameters
        ----------
        r : list[np.ndarray]
            List of coordinate arrays/grids.

        Returns
        -------
        np.ndarray
            The absorbing potential penalty values.

        """
        total_penalty = np.zeros_like(r[0])

        for x, limits in zip(r, self.limits, strict=False):
            lower_penetration = np.maximum(0, limits[0] - x)
            upper_penetration = np.maximum(0, x - limits[1])

            step_penalty = np.exp(self.exponent * (lower_penetration + upper_penetration)) - 1

            total_penalty += step_penalty

        return total_penalty
