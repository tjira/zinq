from ..backend import np
from .grid import Grid
from .hamiltonian import Hamiltonian


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

    def init_gauss(self, grid: Grid, position: np.ndarray, momentum: np.ndarray, gamma: np.ndarray, state: int):
        assert grid.ndim == self.ndim, "INCOMPATIBLE GRID DIMENSION FOR WAVEFUNCTION"
        assert position.shape[0] == self.ndim, "INCOMPATIBLE POSITION VECTOR DIMENSION FOR WAVEFUNCTION"
        assert momentum.shape[0] == self.ndim, "INCOMPATIBLE MOMENTUM VECTOR DIMENSION FOR WAVEFUNCTION"
        assert gamma.shape[0] == self.ndim, "INCOMPATIBLE GAMMA VECTOR DIMENSION FOR WAVEFUNCTION"
        assert state < self.nstate, "STATE INDEX OUT OF BOUNDS FOR WAVEFUNCTION"

        self.data.fill(0)

        exponent = np.zeros_like(grid.position[0], dtype=np.complex128)

        for i in range(self.ndim):
            dr = grid.position[i] - position[i]
            exponent += -0.5 * gamma[i] * dr**2 + 1j * momentum[i] * dr

        self.data[..., state] = np.exp(exponent)
        self.measure = grid.measure
        self.normalize()

    def normalize(self):
        self.data /= np.sqrt(self.norm())

    def diff1(self, grid: Grid, dim: int):
        data_fft = np.fft.fft(self.data, axis=dim)
        d1m = 1j * grid.momentum[dim][..., np.newaxis] * data_fft
        d1x = np.fft.ifft(d1m, axis=dim)
        return Wavefunction.from_data(d1x, self.measure)

    def diff2(self, grid: Grid, dim: int):
        data_fft = np.fft.fft(self.data, axis=dim)
        d2m = -grid.momentum[dim][..., np.newaxis]**2 * data_fft
        d2x = np.fft.ifft(d2m, axis=dim)
        return Wavefunction.from_data(d2x, self.measure)

    def ke(self, grid: Grid, ham: Hamiltonian) -> float:
        return np.sum(self.ke_dens(grid, ham)) * self.measure

    def ke_dens(self, grid: Grid, ham: Hamiltonian) -> np.ndarray:
        ked = np.zeros(self.data.shape[:-1])

        for d in range(self.ndim):
            d2 = self.diff2(grid, d)
            ked += self.overlap_dens(d2).real

        return (-0.5 / ham.mass) * ked

    def momentum(self, grid: Grid) -> np.ndarray:
        return np.array([self.overlap(self.diff1(grid, d)).imag for d in range(self.ndim)])

    def momentum_dens(self, grid: Grid, dim: int) -> np.ndarray:
        return self.overlap_dens(self.diff1(grid, dim)).imag

    def norm(self) -> float:
        return np.vdot(self.data, self.data).real * self.measure

    def norm_dens(self) -> np.ndarray:
        return np.einsum("...i,...i->...", self.data.conj(), self.data).real

    def overlap(self, other) -> complex:
        return np.vdot(self.data, other.data) * self.measure

    def overlap_dens(self, other) -> np.ndarray:
        return np.einsum("...i,...i->...", self.data.conj(), other.data)

    def project_out(self, others):
        for other in others:
            self.data -= self.overlap(other).conj() * other.data

        self.normalize()

    def population(self) -> np.ndarray:
        return np.sum(np.abs(self.data)**2, axis=tuple(range(self.ndim))) * self.measure

    def population_dens(self, state: int) -> np.ndarray:
        return np.abs(self.data[..., state])**2

    def position(self, grid: Grid) -> np.ndarray:
        rho = self.norm_dens()
        return np.array([np.sum(rho * p) for p in grid.position]) * self.measure

    def position_dens(self, grid: Grid, dim: int) -> np.ndarray:
        return self.norm_dens() * grid.position[dim]

    def pe(self, grid: Grid, ham: Hamiltonian, time: float = 0) -> float:
        return np.sum(self.pe_dens(grid, ham, time)) * self.measure

    def pe_dens(self, grid: Grid, ham: Hamiltonian, time: float = 0) -> np.ndarray:
        V = ham.potential.eval_d(grid.position, time=time)
        return np.einsum("...i,ij...,...j->...", self.data.conj(), V, self.data).real

    def to_adiabatic(self, grid: Grid, ham: Hamiltonian) -> 'Wavefunction':
        _, U = np.linalg.eigh( np.moveaxis(ham.potential.eval_d(grid.position), [0, 1], [-2, -1]))
        return Wavefunction.from_data(np.einsum("...ji,...j->...i", U.conj(), self.data), self.measure)
