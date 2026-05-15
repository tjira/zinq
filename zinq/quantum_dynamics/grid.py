import numpy as np


def generate_momentum_grid(limits: np.ndarray, npoint: int):
    assert limits.ndim == 2, f"LIMITS MUST BE 2D ARRAY, GOT {limits.ndim}D"
    assert limits.shape[1] == 2, f"LIMITS MUST HAVE SHAPE (NDIM, 2), GOT {limits.shape}"
    assert npoint > 1, f"NUMBER OF GRID POINTS MUST BE GREATER THAN 1, GOT {npoint}"

    grids = []

    for i in range(limits.shape[0]):
        grids.append(2 * np.pi * np.fft.fftfreq(npoint, d=np.ptp(limits[i]) / (npoint - 1)))

    return np.meshgrid(*grids, indexing="ij")


def generate_position_grid(limits: np.ndarray, npoint: int):
    assert limits.ndim == 2, f"LIMITS MUST BE 2D ARRAY, GOT {limits.ndim}D"
    assert limits.shape[1] == 2, f"LIMITS MUST HAVE SHAPE (NDIM, 2), GOT {limits.shape}"
    assert npoint > 1, f"NUMBER OF GRID POINTS MUST BE GREATER THAN 1, GOT {npoint}"

    grids = []

    for i in range(limits.shape[0]):
        grids.append(np.linspace(limits[i, 0], limits[i, 1], npoint))

    return np.meshgrid(*grids, indexing="ij")
