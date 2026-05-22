from functools import cached_property

import numpy as np

from .grid import Grid
from .hamiltonian import Hamiltonian
from .initial_conditions import InitialConditions


class Wavefunction:
    def __init__(self, ic: InitialConditions, grid: Grid, H: Hamiltonian):
        self.data = np.zeros((*[grid.npoint] * grid.ndim, H.pot.nstate), dtype=np.complex128)

        exponent = np.zeros_like(grid.pos[0], dtype=np.complex128)

        for pos, r, k, g in zip(grid.pos, ic.pos, ic.mom, ic.gamma):
            exponent += -0.5 * g * (dr := pos - r)**2 + 1j * k * dr

        self.data[..., ic.state], self.absorbed = np.exp(exponent), np.zeros(H.pot.nstate)

        self.normalize(grid)

        if ic.adia:
            self.data = self.to_dia(H).data

    @classmethod
    def from_data(cls, data: np.ndarray, absorbed: np.ndarray | None = None):
        (wfn := cls.__new__(cls)).__dict__.update({"data": data, "absorbed": np.zeros(data.shape[-1])})

        if absorbed is not None:
            wfn.absorbed = absorbed

        return wfn

    @cached_property
    def ndim(self) -> int:
        return len(self.data.shape) - 1

    @cached_property
    def nstate(self) -> int:
        return self.data.shape[-1]

    @property
    def density(self) -> np.ndarray:
        return np.abs(self.data)**2

    @property
    def rho(self) -> np.ndarray:
        return np.sum(self.density, axis=-1)

    def copy(self) -> "Wavefunction":
        return Wavefunction.from_data(self.data.copy(), self.absorbed.copy())

    def ke(self, grid: Grid, H: Hamiltonian) -> float:
        return np.sum(self.to_kspace().rho * H.T) * grid.measure

    def mom(self, grid: Grid) -> np.ndarray:
        rho_k = self.to_kspace().rho

        return np.array([np.sum(rho_k * k) for k in grid.mom]) * grid.measure

    def norm(self, grid: Grid) -> float:
        return np.real(np.vdot(self.data, self.data)).item() * grid.measure

    def normalize(self, grid: Grid):
        self.data /= np.sqrt(norm) if (norm := self.norm(grid)) > 1e-14 else 1

    def overlap(self, other: "Wavefunction", grid: Grid) -> complex:
        return np.vdot(self.data, other.data).item() * grid.measure

    def pe(self, grid: Grid, H: Hamiltonian) -> float:
        pe_dens = np.einsum("...i,...ij,...j->...", np.conj(self.data), H.V, self.data)

        return np.real(np.sum(pe_dens)) * grid.measure

    def pop(self, grid: Grid) -> np.ndarray:
        return np.sum(self.density, axis=tuple(range(grid.ndim))) * grid.measure

    def pos(self, grid: Grid) -> np.ndarray:
        rho = self.rho

        return np.array([np.sum(rho * r) for r in grid.pos]) * grid.measure

    def project_out(self, others: list["Wavefunction"], grid: Grid):
        for other in others:
            self.data -= np.conj(self.overlap(other, grid)) * other.data

        self.normalize(grid)

    def to_adia(self, H: Hamiltonian) -> "Wavefunction":
        return Wavefunction.from_data(np.einsum("...ji,...j->...i", np.conj(H.U), self.data), self.absorbed)

    def to_dia(self, H: Hamiltonian) -> "Wavefunction":
        return Wavefunction.from_data(np.einsum("...ij,...j->...i", H.U, self.data), self.absorbed)

    def to_kspace(self) -> "Wavefunction":
        return Wavefunction.from_data(np.fft.fftn(self.data, axes=range(self.ndim), norm="ortho"), self.absorbed)
