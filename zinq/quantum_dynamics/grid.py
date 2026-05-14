import numpy as np

def generateMomentumGrid(limits: np.ndarray, npoint: int):
    assert limits.shape[1] == 2

    grids = []

    for i in range(limits.shape[0]):
        grids.append(2 * np.pi * np.fft.fftfreq(npoint, d=np.ptp(limits[i]) / (npoint - 1)))

    return np.meshgrid(*grids, indexing="ij")

def generatePositionGrid(limits: np.ndarray, npoint: int):
    assert limits.shape[1] == 2

    grids = []

    for i in range(limits.shape[0]):
        grids.append(np.linspace(limits[i, 0], limits[i, 1], npoint))

    return np.meshgrid(*grids, indexing="ij")
