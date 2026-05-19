import numpy as np
from typing_extensions import TypedDict


class Grid:
    Config = TypedDict("Config", {
        "limits": list[tuple[float, float]],
        "npoint": int,
    })

    position: list[np.ndarray]
    momentum: list[np.ndarray]

    @property
    def ndim(self) -> int:
        return len(self.position)

    @property
    def npoint(self) -> int:
        return self.position[0].shape[0]

    @property
    def limits(self) -> np.ndarray:
        return np.array([[p.min(), p.max() + np.ptp(p) / (self.npoint - 1)] for p in self.position])

    @property
    def measure(self) -> float:
        return np.prod(np.ptp(self.limits, axis=1) / self.npoint)

    def __init__(self, limits: np.ndarray, npoint: int):
        self.position = self._position_grid(np.array(limits), npoint)
        self.momentum = self._momentum_grid(np.array(limits), npoint)

    def _momentum_grid(self, limits: np.ndarray, npoint: int) -> list[np.ndarray]:
        grids = []

        for i in range(limits.shape[0]):
            grids.append(2 * np.pi * np.fft.fftfreq(npoint, d=np.ptp(limits[i]) / npoint))

        return list(np.meshgrid(*grids, indexing="ij"))

    def _position_grid(self, limits: np.ndarray, npoint: int) -> list[np.ndarray]:
        grids = []

        for i in range(limits.shape[0]):
            grids.append(np.linspace(limits[i, 0], limits[i, 1], npoint, endpoint=False))

        return list(np.meshgrid(*grids, indexing="ij"))
