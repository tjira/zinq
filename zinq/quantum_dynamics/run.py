import datetime
import time

import numpy as np

from ..potential import Potential
from . import strang_split
from .grid import generateMomentumGrid, generatePositionGrid
from .options import Options
from .strang_split import StrangSplit
from .wavefunction import Wavefunction


def run(options_dict: dict):
    opt = Options(**options_dict)

    potential = opt.potential.create()

    wfn = Wavefunction(potential.ndim, potential.nstate, opt.grid.npoint)

    position_grid = generatePositionGrid(np.array(opt.grid.limits), opt.grid.npoint)
    momentum_grid = generateMomentumGrid(np.array(opt.grid.limits), opt.grid.npoint)

    wfn.initializeGaussian(
        position_grid=position_grid,
        position=np.array(opt.initial_conditions.position),
        momentum=np.array(opt.initial_conditions.momentum),
        gamma=np.array(opt.initial_conditions.gamma),
        state=opt.initial_conditions.state
    )

    propagator = StrangSplit(
        position_grid=position_grid,
        momentum_grid=momentum_grid,
        potential=potential,
        mass=opt.mass,
        dt=opt.time_step,
        imaginary=opt.imaginary,
    )

    with np.printoptions(formatter={"float": "{:10.4f}".format}, suppress=True):
        print(f"\nINITIAL GAMMA: {np.array(opt.initial_conditions.gamma, dtype=float)}\n")

    print(
        f"{'ITER':>5} {'KIN (Eh)':>12} {'POT (Eh)':>12} {'TOT (Eh)':>12} "
        f"{'POS (a0)':>{(11 * wfn.ndim + 1)}} {'MOM (hb/a0)':>{(11 * wfn.ndim + 1)}} "
        f"{'POPULATION':>{(11 * wfn.nstate + 1)}} {'NORM':>9} TIME"
    )

    for i in range(opt.iterations + 1):
        start_time = time.time()

        if i: propagator.step(wfn)

        log_iteration = i == 0 or i == opt.iterations or (i % opt.log_interval == 0)

        if not log_iteration: continue

        kinetic_energy = wfn.kineticEnergy(momentum_grid, opt.mass)
        potential_energy = wfn.potentialEnergy(position_grid, potential)
        total_energy = kinetic_energy + potential_energy
        position = wfn.position(position_grid)
        momentum = wfn.momentum(momentum_grid)
        population = wfn.population()
        norm = wfn.norm()

        duration = datetime.timedelta(seconds=time.time() - start_time)

        with np.printoptions(formatter={"float": "{:10.4f}".format}, suppress=True):
            print(
                f"{i:5d} {kinetic_energy:12.6f} {potential_energy:12.6f} {total_energy:12.6f} "
                f"{position} {momentum} {population} {norm:1.3e} {duration}"
            )
