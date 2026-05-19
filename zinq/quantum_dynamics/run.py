import numpy as np

from .initial_conditions import InitialConditions
from ..zig import load_library
from .grid import Grid
from .wavefunction import Wavefunction
from .hamiltonian import Hamiltonian
from .options import Options


def run(opt: Options):
    grid, wfn, ham = _init(opt)
    print(f"ZIG ADD(10, 32) = {load_library('native').add(10, 32)}")
    print(wfn.data)

def _init(opt: Options) -> tuple[Grid, Wavefunction, Hamiltonian]:
    grid = Grid(
        limits=np.array(opt.grid["limits"]),
        npoint=opt.grid["npoint"],
    )
    wfn = Wavefunction(
        ic=InitialConditions(**opt.initial_conditions),
        grid=grid,
        nstate=opt.potential.nstate,
    )
    ham = Hamiltonian(
        pot=opt.potential,
        cap=None,
        mass=opt.mass,
    )

    return grid, wfn, ham
