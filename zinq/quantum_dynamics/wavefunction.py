import numpy as np

from ..potential import Potential


class Wavefunction:
    data: np.ndarray
    measure: float

    def __init__(self, ndim: int, nstate: int, npoint: int):
        self.data = np.empty((*[npoint] * ndim, nstate), dtype=np.complex128)

    @classmethod
    def from_data(cls, data: np.ndarray, measure: float):
        wfn = cls.__new__(cls)
        
        wfn.data = data
        wfn.measure = measure
        
        return wfn

    @property
    def ndim(self) -> int:
        return self.data.ndim - 1

    @property
    def nstate(self) -> int:
        return self.data.shape[-1]

    @property
    def npoint(self) -> int:
        return self.data.shape[0]

    def initialize_gaussian(self, /,
        position_grid: list[np.ndarray],
        position: np.ndarray,
        momentum: np.ndarray,
        gamma: np.ndarray,
        state: int
    ):
        assert len(position_grid) == self.ndim, (
            f"POSITION GRID DIMENSION DOES NOT MATCH WAVEFUNCTION DIMENSION"
        )
        assert position.shape[0] == self.ndim, (
            f"POSITION VECTOR DIMENSION DOES NOT MATCH WAVEFUNCTION DIMENSION"
        )
        assert momentum.shape[0] == self.ndim, (
            f"MOMENTUM VECTOR DIMENSION DOES NOT MATCH WAVEFUNCTION DIMENSION"
        )
        assert gamma.shape[0] == self.ndim, (
            f"GAMMA VECTOR DIMENSION DOES NOT MATCH WAVEFUNCTION DIMENSION"
        )
        assert state < self.nstate, (
            f"INITIAL STATE ({state}) IS OUT OF BOUNDS FOR {self.nstate} STATES"
        )

        self.data.fill(0)

        exponent = np.zeros_like(position_grid[0], dtype=np.complex128)

        for i in range(self.ndim):
            dr = position_grid[i] - position[i]
            exponent += -0.5 * gamma[i] * dr**2 + 1j * momentum[i] * dr

        self.data[..., state] = np.exp(exponent)
        self.measure = np.prod([np.ptp(axis) / (self.npoint - 1) for axis in position_grid])
        self.normalize()

    def normalize(self):
        self.data /= np.sqrt(self.norm())

    def first_derivative(self, momentum_grid: list[np.ndarray], dim: int):
        data_fft = np.fft.fft(self.data, axis=dim)
        first_derivative_momentum = 1j * momentum_grid[dim][..., np.newaxis] * data_fft
        first_derivative_position = np.fft.ifft(first_derivative_momentum, axis=dim)
        return Wavefunction.from_data(first_derivative_position, self.measure)

    def second_derivative(self, momentum_grid: list[np.ndarray], dim: int):
        data_fft = np.fft.fft(self.data, axis=dim)
        second_derivative_momentum = -momentum_grid[dim][..., np.newaxis]**2 * data_fft
        second_derivative_position = np.fft.ifft(second_derivative_momentum, axis=dim)
        return Wavefunction.from_data(second_derivative_position, self.measure)

    def kinetic_energy(self, momentum_grid: list[np.ndarray], mass: float) -> float:
        return np.sum(self.kinetic_energy_density(momentum_grid, mass)) * self.measure

    def kinetic_energy_density(self, momentum_grid: list[np.ndarray], mass: float) -> np.ndarray:
        kinetic_energy_density = np.zeros(self.data.shape[:-1])

        for d in range(self.ndim):
            second_derivative = self.second_derivative(momentum_grid, d)
            kinetic_energy_density += self.overlap_density(second_derivative).real

        return (-0.5 / mass) * kinetic_energy_density

    def momentum(self, momentum_grid: list[np.ndarray]) -> np.ndarray:
        momentum = np.zeros(self.ndim)

        for d in range(self.ndim):
            momentum[d] = np.sum(self.momentum_density(momentum_grid, d)) * self.measure

        return momentum

    def momentum_density(self, momentum_grid: list[np.ndarray], dim: int) -> np.ndarray:
        return self.overlap_density(self.first_derivative(momentum_grid, dim)).imag

    def norm(self) -> float:
        return np.sum(self.norm_density()) * self.measure

    def norm_density(self) -> np.ndarray:
        return np.sum(np.abs(self.data)**2, axis=-1)

    def overlap(self, other) -> complex:
        return np.sum(self.overlap_density(other)) * self.measure

    def overlap_density(self, other) -> np.ndarray:
        return np.einsum("...i,...i->...", np.conj(self.data), other.data)

    def project_out(self, others):
        for other in others:
            self.data -= np.conj(self.overlap(other)) * other.data

        self.normalize()

    def population(self) -> np.ndarray:
        population = np.zeros(self.nstate)

        for s in range(self.nstate):
            population[s] = np.sum(self.population_density(s)) * self.measure

        return population

    def population_density(self, state: int) -> np.ndarray:
        return np.abs(self.data[..., state])**2

    def position(self, position_grid: list[np.ndarray]) -> np.ndarray:
        position = np.zeros(self.ndim)

        for d in range(self.ndim):
            position[d] = np.sum(self.position_density(position_grid, d)) * self.measure

        return position

    def position_density(self, position_grid: list[np.ndarray], dim: int) -> np.ndarray:
        return np.sum(np.abs(self.data)**2, axis=-1) * position_grid[dim]

    def potential_energy(self,
        position_grid: list[np.ndarray],
        potential: Potential,
        time: float = 0
    ) -> float:
        potential_energy_density = self.potential_energy_density(position_grid, potential, time)
        return np.sum(potential_energy_density) * self.measure

    def potential_energy_density(self,
        position_grid: list[np.ndarray],
        potential: Potential,
        time: float = 0
    ) -> np.ndarray:
        V = potential.evaluate_diabatic(position_grid, time=time)
        return np.einsum("...i,ij...,...j->...", np.conj(self.data), V, self.data).real
