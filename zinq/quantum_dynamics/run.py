import datetime, time

import numpy as np

from .grid import generatePositionGrid
from .grid import generateMomentumGrid
from .split_operator import SplitOperator
from .wavefunction import Wavefunction

from ..potential import TullyFirst, Harmonic

def run(options: dict):
    potential = Harmonic(k=np.array(options["potential"]["harmonic"]["k"]))
    # pot = TullyFirst()

    wfn = Wavefunction(potential.ndim, potential.nstate, options["grid"]["npoint"])

    position_grid = generatePositionGrid(np.array(options["grid"]["limits"]), options["grid"]["npoint"])
    momentum_grid = generateMomentumGrid(np.array(options["grid"]["limits"]), options["grid"]["npoint"])

    wfn.initializeGaussian(
        position_grid,
        np.array(options["initial_conditions"]["position"]),
        np.array(options["initial_conditions"]["momentum"]),
        np.array(options["initial_conditions"]["gamma"]),
        np.array(options["initial_conditions"]["state"])
    )

    propagator = SplitOperator(position_grid, momentum_grid, potential, options["time_step"], options["mass"], True)

    with np.printoptions(formatter={"float": "{:10.4f}".format}, suppress=True):
        print(f"\nINITIAL GAMMA: {np.array(options['initial_conditions']['gamma'], dtype=float)}\n")

    print(
        f"{'ITER':>5} {'KIN (Eh)':>12} {'POT (Eh)':>12} {'TOT (Eh)':>12} "
        f"{'POS (a0)':>{(11 * wfn.ndim + 1)}} {'MOM (hb/a0)':>{(11 * wfn.ndim + 1)}} "
        f"{'POPULATION':>{(11 * wfn.nstate + 1)}} {'NORM':>9} TIME"
    )

    for i in range(options["iterations"] + 1):
        start_time = time.time()

        if i: propagator.step(wfn)

        log_iteration = i == 0 or i == options["iterations"] or (i % options.get("log_interval", 1) == 0)

        if log_iteration:
            kinetic_energy = wfn.kineticEnergy(momentum_grid, options["mass"])
            potential_energy = wfn.potentialEnergy(position_grid, potential)
            total_energy = kinetic_energy + potential_energy
            position = wfn.position(position_grid)
            momentum = wfn.momentum(momentum_grid)
            population = wfn.population()
            norm = wfn.norm()

        if not log_iteration: continue

        duration = datetime.timedelta(seconds=time.time() - start_time)

        with np.printoptions(formatter={"float": "{:10.4f}".format}, suppress=True):
            print(
                f"{i:5d} {kinetic_energy:12.6f} {potential_energy:12.6f} {total_energy:12.6f} "
                f"{position} {momentum} {population} {norm:1.3e} {duration}"
            )
