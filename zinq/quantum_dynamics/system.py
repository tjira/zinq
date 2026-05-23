from dataclasses import dataclass

import numpy as np

from .grid import Grid
from .hamiltonian import Hamiltonian
from .absorber import Absorber
from .initial_conditions import InitialConditions
from .wavefunction import Wavefunction
from .options import Options


@dataclass
class System:
    ham: Hamiltonian
    grid: Grid
    wfn: Wavefunction

    @classmethod
    def from_options(cls, opt: Options):
        grid, absorber = Grid.from_options(opt.grid), None

        if opt.hamiltonian.absorber:
            absorber = Absorber.from_options(opt.hamiltonian.absorber)

        ham = Hamiltonian(
            grid=grid,
            pot=opt.hamiltonian.potential,
            m=opt.hamiltonian.mass,
            absorber=absorber,
        )

        wfn = Wavefunction(
            ic=InitialConditions.from_options(opt.initial_conditions),
            grid=grid,
            ham=ham,
        )

        return cls(ham=ham, grid=grid, wfn=wfn)

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
