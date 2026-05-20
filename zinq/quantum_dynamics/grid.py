from functools import cached_property

import numpy as np

from .options import GridConfig


class Grid:
    pos: list[np.ndarray]
    mom: list[np.ndarray]

    def __init__(self, limits: np.ndarray, npoint: int):
        self.pos = self._pos_grid(np.array(limits), npoint)
        self.mom = self._mom_grid(np.array(limits), npoint)

    @classmethod
    def from_options(cls, opt: GridConfig):
        return cls(limits=np.array(opt.limits), npoint=opt.npoint)

    @cached_property
    def ndim(self) -> int:
        return len(self.pos)

    @cached_property
    def npoint(self) -> int:
        return self.pos[0].shape[0]

    @cached_property
    def limits(self) -> np.ndarray:
        return np.array([[p.min(), p.max()] for p in self.pos])

    @cached_property
    def measure(self) -> float:
        return np.prod(np.ptp(self.limits, axis=1) / (self.npoint - 1))

    def _mom_grid(self, limits: np.ndarray, npoint: int) -> list[np.ndarray]:
        grids = []

        for i in range(limits.shape[0]):
            grids.append(2 * np.pi * np.fft.fftfreq(npoint, d=np.ptp(limits[i]) / (npoint - 1)))

        return list(np.meshgrid(*grids, indexing="ij"))

    def _pos_grid(self, limits: np.ndarray, npoint: int) -> list[np.ndarray]:
        grids = []

        for i in range(limits.shape[0]):
            grids.append(np.linspace(limits[i, 0], limits[i, 1], npoint, endpoint=True))

        return list(np.meshgrid(*grids, indexing="ij"))
