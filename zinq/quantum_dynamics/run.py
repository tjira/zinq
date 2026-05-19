from itertools import count
import time

import numpy as np

from ..potential import Potential
from ..zig import load_library
from .grid import Grid
from .hamiltonian import Hamiltonian
from .initial_conditions import InitialConditions
from .options import Options
from .strang_split import StrangSplit
from .wavefunction import Wavefunction


def run(opt: Options):
    grid, wfn, H, pot, prop, start_time = *_init(opt), time.time()

    for i in range(opt.iterations + 1) if opt.iterations else count():
        current_time = i * prop.dt

        if i: prop.step(wfn, grid, H, pot, current_time - prop.dt)


        norm = wfn.norm(grid)
        pop = wfn.pop(grid)
        pos = wfn.pos(grid)
        mom = wfn.mom(grid)
        pe = wfn.pe(grid, H)
        ke = wfn.ke(grid, H)
        # print(norm)
        print(pop)
        # print(pos)
        # print(mom)
        # print(pe)
        # print(ke)
    # print(f"ZIG ADD(10, 32) = {load_library('native').add(10, 32)}")

def _init(opt: Options) -> tuple[Grid, Wavefunction, Hamiltonian, Potential, StrangSplit]:
    grid = Grid(
        limits=np.array(opt.grid.limits),
        npoint=opt.grid.npoint,
    )
    H = Hamiltonian(
        grid=grid,
        pot=opt.hamiltonian.potential,
        m=opt.hamiltonian.mass,
    )
    ic = InitialConditions(
        pos=np.array(opt.initial_conditions.position),
        mom=np.array(opt.initial_conditions.momentum),
        gamma=np.array(opt.initial_conditions.gamma),
        state=opt.initial_conditions.state,
        adia=opt.initial_conditions.adiabatic,
    )
    wfn = Wavefunction(
        ic=ic,
        grid=grid,
        nstate=opt.hamiltonian.potential.nstate,
    )
    prop = StrangSplit(
        grid=grid,
        H=H,
        pot=opt.hamiltonian.potential,
        dt=opt.time_step,
        imaginary=opt.imaginary is not None,
    )
    return grid, wfn, H, opt.hamiltonian.potential, prop
