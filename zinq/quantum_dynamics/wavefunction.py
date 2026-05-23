from functools import cached_property

import numpy as np

from .grid import Grid
from .hamiltonian import Hamiltonian
from .initial_conditions import InitialConditions
from .options import InitialConditionsConfig


class Wavefunction:
    def __init__(self, *, ic: InitialConditions, grid: Grid, ham: Hamiltonian):
        self.data = np.zeros((*[grid.npoint] * grid.ndim, ham.pot.nstate), dtype=np.complex128)

        exponent = np.zeros_like(grid.pos[0], dtype=np.complex128)

        for pos, r, k, g in zip(grid.pos, ic.pos, ic.mom, ic.gamma):
            exponent += -0.5 * g * (dr := pos - r)**2 + 1j * k * dr

        self.data[..., ic.state] = np.exp(exponent)

        self.normalize(grid)

        if ic.adia:
            self.data = self.to_dia(ham).data

    @classmethod
    def from_data(cls, data: np.ndarray):
        return (wfn := cls.__new__(cls)).__dict__.update({"data": data}) or wfn

    @classmethod
    def from_options(cls, opt: InitialConditionsConfig, grid: Grid, ham: Hamiltonian):
        return cls(ic=InitialConditions.from_options(opt), grid=grid, ham=ham)

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
        return Wavefunction.from_data(self.data.copy())

    def ke(self, grid: Grid, ham: Hamiltonian) -> float:
        return np.sum(self.to_kspace().rho * ham.T) * grid.measure

    def mom(self, grid: Grid) -> np.ndarray:
        rho_k = self.to_kspace().rho

        return np.array([np.sum(rho_k * k) for k in grid.mom]) * grid.measure

    def norm(self, grid: Grid) -> float:
        return np.real(np.vdot(self.data, self.data)).item() * grid.measure

    def normalize(self, grid: Grid):
        self.data /= np.sqrt(norm) if (norm := self.norm(grid)) > 1e-14 else 1

    def overlap(self, other: "Wavefunction", grid: Grid) -> complex:
        return np.vdot(self.data, other.data).item() * grid.measure

    def pe(self, grid: Grid, ham: Hamiltonian) -> float:
        pe_dens = np.einsum("...i,...ij,...j->...", np.conj(self.data), ham.V, self.data)

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

    def to_adia(self, ham: Hamiltonian) -> "Wavefunction":
        return Wavefunction.from_data(np.einsum("...ji,...j->...i", np.conj(ham.U), self.data))

    def to_dia(self, ham: Hamiltonian) -> "Wavefunction":
        return Wavefunction.from_data(np.einsum("...ij,...j->...i", ham.U, self.data))

    def to_kspace(self) -> "Wavefunction":
        return Wavefunction.from_data(np.fft.fftn(self.data, axes=range(self.ndim), norm="ortho"))
