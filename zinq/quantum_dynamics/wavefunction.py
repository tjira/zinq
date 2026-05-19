import numpy as np

from .grid import Grid
from .hamiltonian import Hamiltonian
from .initial_conditions import InitialConditions


class Wavefunction:
    data: np.ndarray

    def __init__(self, ic: InitialConditions, grid: Grid, nstate: int):
        self.data = np.zeros((*[grid.npoint] * grid.ndim, nstate), dtype=np.complex128)
        exponent = np.zeros_like(grid.pos[0], dtype=np.complex128)

        for pos, r, k, g in zip(grid.pos, ic.pos, ic.mom, ic.gamma):
            exponent += -0.5 * g * (dr := pos - r)**2 + 1j * k * dr

        self.data[..., ic.state] = np.exp(exponent)

        self.normalize(grid)

    def ke(self, grid: Grid, H: Hamiltonian) -> float:
        data_k = np.fft.fftn(self.data, axes=range(grid.ndim), norm="ortho")
        ke_dens = np.real(np.einsum("...i,...,...i->...", np.conj(data_k), H.T, data_k))

        return np.sum(ke_dens) * grid.measure

    def mom(self, grid: Grid) -> np.ndarray:
        data_k = np.fft.fftn(self.data, axes=range(grid.ndim), norm="ortho")
        rho_k = np.sum(np.abs(data_k)**2, axis=-1)

        return np.array([np.sum(rho_k * k) for k in grid.mom]) * grid.measure

    def norm(self, grid: Grid) -> float:
        return np.real(np.vdot(self.data, self.data)) * grid.measure

    def normalize(self, grid: Grid):
        self.data /= np.sqrt(self.norm(grid))

    def pe(self, grid: Grid, H: Hamiltonian) -> float:
        pe_dens = np.real(np.einsum("...i,...ij,...j->...", np.conj(self.data), H.V, self.data))

        return np.sum(pe_dens) * grid.measure

    def pop(self, grid: Grid) -> np.ndarray:
        return np.sum(np.abs(self.data)**2, axis=tuple(range(grid.ndim))) * grid.measure

    def pos(self, grid: Grid) -> np.ndarray:
        rho = np.sum(np.abs(self.data)**2, axis=-1)

        return np.array([np.sum(rho * r) for r in grid.pos]) * grid.measure
