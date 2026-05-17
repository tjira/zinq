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
        norm = self.norm()
        if norm > 1e-14:
            self.data /= np.sqrt(norm)

    def ke(self, grid: Grid, H: Hamiltonian) -> float:
        data_fft = np.fft.fftn(self.data, axes=range(self.ndim))
        k_sq = sum((k**2 for k in grid.momentum), np.zeros_like(grid.momentum[0]))
        factor = (0.5 / H.mass) * self.measure / np.prod(self.data.shape[:-1])

        return float(np.sum(np.abs(data_fft)**2 * k_sq[..., np.newaxis]) * factor)

    def momentum(self, grid: Grid) -> np.ndarray:
        data_fft = np.fft.fftn(self.data, axes=range(self.ndim))
        abs_sq = np.abs(data_fft)**2
        factor = self.measure / np.prod(self.data.shape[:-1])

        return np.array([np.sum(abs_sq * k[..., np.newaxis]) * factor for k in grid.momentum])

    def norm(self) -> float:
        return float(np.real(np.vdot(self.data, self.data)) * self.measure)

    def overlap(self, other) -> complex:
        return complex(np.vdot(self.data, other.data) * self.measure)

    def project_out(self, others):
        for other in others:
            self.data -= np.conj(self.overlap(other)) * other.data

        self.normalize()

    def population(self, pop_decay: np.ndarray) -> np.ndarray:
        return np.sum(np.abs(self.data)**2, axis=tuple(range(self.ndim))) * self.measure + pop_decay

    def position(self, grid: Grid) -> np.ndarray:
        rho = np.real(np.einsum("...i,...i->...", np.conj(self.data), self.data))

        return np.array([np.sum(rho * p) for p in grid.position]) * self.measure

    def pe(self, grid: Grid, H: Hamiltonian, time: float = 0) -> float:
        V = H.pot.eval_d(grid.position, time=time)
        pe_dens = np.real(np.einsum("...i,ij...,...j->...", np.conj(self.data), V, self.data))

        return float(np.sum(pe_dens) * self.measure)

    def to_adiabatic(self, grid: Grid, H: Hamiltonian, time: float = 0) -> 'Wavefunction':
        _, U = np.linalg.eigh(np.moveaxis(H.pot.eval_d(grid.position, time), [0, 1], [-2, -1]))
        return Wavefunction.from_data(np.einsum("...ji,...j->...i", np.conj(U), self.data), self.measure)

    def to_diabatic(self, grid: Grid, H: Hamiltonian, time: float = 0) -> 'Wavefunction':
        _, U = np.linalg.eigh(np.moveaxis(H.pot.eval_d(grid.position, time), [0, 1], [-2, -1]))
        return Wavefunction.from_data(np.einsum("...ij,...j->...i", U, self.data), self.measure)
