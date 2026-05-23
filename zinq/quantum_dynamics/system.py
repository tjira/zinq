from dataclasses import dataclass

import numpy as np

from .grid import Grid
from .hamiltonian import Hamiltonian
from .wavefunction import Wavefunction
from .options import Options


@dataclass
class System:
    ham: Hamiltonian
    grid: Grid
    wfn: Wavefunction

    @classmethod
    def from_options(cls, opt: Options):
        return cls(
            grid=(grid := Grid.from_options(opt.grid)),
            ham=(ham := Hamiltonian.from_options(opt.hamiltonian, grid)),
            wfn=Wavefunction.from_options(opt.initial_conditions, grid, ham),
        )

    @property
    def ke(self) -> float:
        return self.wfn.ke(self.grid, self.ham)

    @property
    def mom(self) -> np.ndarray:
        return self.wfn.mom(self.grid)

    @property
    def norm(self) -> float:
        return self.wfn.norm(self.grid)

    @property
    def pe(self) -> float:
        return self.wfn.pe(self.grid, self.ham)

    @property
    def pop(self) -> np.ndarray:
        return self.wfn.pop(self.grid)

    @property
    def pos(self) -> np.ndarray:
        return self.wfn.pos(self.grid)

    def normalize(self):
        self.wfn.normalize(self.grid)

    def update_V(self, time: float):
        self.ham.update_V(self.grid, time)
