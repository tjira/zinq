from pydantic import BaseModel

from ..backend import np


class GridOptions(BaseModel):
    limits: list[list[float]]
    npoint: int


class Grid:
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
        return np.array([[p.min(), p.max()] for p in self.position])

    @property
    def measure(self) -> float:
        return np.prod(np.ptp(self.limits, axis=1) / (self.npoint - 1))

    def __init__(self, limits: np.ndarray, npoint: int):
        assert limits.ndim == 2, f"LIMITS MUST BE 2D ARRAY, GOT {limits.ndim}D"
        assert limits.shape[1] == 2, f"LIMITS MUST HAVE SHAPE (NDIM, 2), GOT {limits.shape}"
        assert npoint > 1, f"NUMBER OF GRID POINTS MUST BE GREATER THAN 1, GOT {npoint}"

        self.position = self._pos_grid(np.array(limits), npoint)
        self.momentum = self._mom_grid(np.array(limits), npoint)

    def _mom_grid(self, limits: np.ndarray, npoint: int) -> list[np.ndarray]:
        grids = []
        for i in range(limits.shape[0]):
            grids.append(2 * np.pi * np.fft.fftfreq(npoint, d=np.ptp(limits[i]) / (npoint - 1)))
        return list(np.meshgrid(*grids, indexing="ij"))

    def _pos_grid(self, limits: np.ndarray, npoint: int) -> list[np.ndarray]:
        grids = []
        for i in range(limits.shape[0]):
            grids.append(np.linspace(limits[i, 0], limits[i, 1], npoint))
        return list(np.meshgrid(*grids, indexing="ij"))
