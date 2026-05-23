"""Module representing the spatial grid for quantum dynamics simulations."""

from functools import cached_property

import numpy as np

from .options import GridConfig


class Grid:
    """
    Representation of a multidimensional spatial grid and its momentum-space counterpart.

    Attributes
    ----------
    pos : list[np.ndarray]
        List of coordinate meshes along each dimension.
    mom : list[np.ndarray]
        List of momentum meshes along each dimension.

    """

    pos: list[np.ndarray]
    mom: list[np.ndarray]

    def __init__(self, limits: np.ndarray, npoint: int) -> None:
        """
        Initialize the coordinate and momentum grids.

        Parameters
        ----------
        limits : np.ndarray
            Coordinate boundaries for each dimension. Shape (ndim, 2).
        npoint : int
            Number of grid points per dimension.

        """
        self.pos = self._pos_grid(np.array(limits), npoint)
        self.mom = self._mom_grid(np.array(limits), npoint)

    @classmethod
    def from_options(cls, opt: GridConfig) -> "Grid":
        """
        Create a Grid instance from a GridConfig configuration object.

        Parameters
        ----------
        opt : GridConfig
            The grid configuration options.

        Returns
        -------
        Grid
            The initialized Grid instance.

        """
        return cls(limits=np.array(opt.limits), npoint=opt.npoint)

    @cached_property
    def ndim(self) -> int:
        """
        Get the dimensionality of the grid.

        Returns
        -------
        int
            Number of spatial dimensions.

        """
        return len(self.pos)

    @cached_property
    def npoint(self) -> int:
        """
        Get the number of grid points per dimension.

        Returns
        -------
        int
            Number of points per dimension.

        """
        return self.pos[0].shape[0]

    @cached_property
    def limits(self) -> np.ndarray:
        """
        Get the coordinate limits of the grid.

        Returns
        -------
        np.ndarray
            Limits of the grid coordinates for each dimension.

        """
        return np.array([[p.min(), p.max() + np.ptp(p) / (self.npoint - 1)] for p in self.pos])

    @cached_property
    def measure(self) -> float:
        """
        Get the integration measure (volume element) of the grid.

        Returns
        -------
        float
            The spatial volume element.

        """
        return np.prod(np.ptp(self.limits, axis=1) / self.npoint)

    def _mom_grid(self, limits: np.ndarray, npoint: int) -> list[np.ndarray]:
        """Generate momentum grid meshes."""
        grids = [
            2 * np.pi * np.fft.fftfreq(npoint, d=np.ptp(limits[i]) / npoint)
            for i in range(limits.shape[0])
        ]
        return list(np.meshgrid(*grids, indexing="ij"))

    def _pos_grid(self, limits: np.ndarray, npoint: int) -> list[np.ndarray]:
        """Generate spatial grid meshes."""
        grids = [
            np.linspace(limits[i, 0], limits[i, 1], npoint, endpoint=False)
            for i in range(limits.shape[0])
        ]
        return list(np.meshgrid(*grids, indexing="ij"))
