import datetime, time

import numpy as np

from .grid import generatePositionGrid
from .grid import generateMomentumGrid
from .split_operator import SplitOperator
from .wavefunction import Wavefunction

from ..potential import TullyFirst

def run(options: dict):
    pot = TullyFirst();

    wfn = Wavefunction(pot.ndim, pot.nstate, options["grid"]["npoint"])
    propagator = SplitOperator(options["time_step"], options["mass"], False)

    position_grid = generatePositionGrid(np.array(options["grid"]["limits"]), options["grid"]["npoint"])
    momentum_grid = generateMomentumGrid(np.array(options["grid"]["limits"]), options["grid"]["npoint"])

    wfn.initializeGaussian(
        position_grid,
        np.array(options["initial_conditions"]["position"]),
        np.array(options["initial_conditions"]["momentum"]),
        np.array(options["initial_conditions"]["gamma"]),
        np.array(options["initial_conditions"]["state"])
    )

    with np.printoptions(formatter={"float": "{:10.4f}".format}, suppress=True):
        print(f"\nINITIAL GAMMA: {np.array(options['initial_conditions']['gamma'], dtype=float)}\n")

    print(
        f"{'ITER':>5} {'KIN (Eh)':>12} {'POT (Eh)':>12} {'TOT (Eh)':>12} "
        f"{'POS (a0)':>{(11 * wfn.ndim + 1)}} {'MOM (hb/a0)':>{(11 * wfn.ndim + 1)}} "
        f"{'POPULATION':>{(11 * wfn.nstate + 1)}} {'NORM':>9} TIME"
    )

    for i in range(options["iterations"] + 1):
        start_time = time.time()

        if i: propagator.step(wfn, position_grid, momentum_grid, pot)

        kinetic_energy = wfn.kineticEnergy(momentum_grid, options["mass"])
        potential_energy = wfn.potentialEnergy(position_grid, pot)
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

    # import matplotlib.pyplot as plt
    # plt.plot(grid[0], wfn.data.real[..., 1], label="state 0")
    # plt.plot(grid[0], wfn.data.imag[..., 1], label="state 1")
    # plt.show()
