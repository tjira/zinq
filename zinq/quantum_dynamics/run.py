from dataclasses import dataclass
import datetime
import time

import numpy as np

from ..potential import Potential
from . import strang_split
from .grid import generateMomentumGrid, generatePositionGrid
from .options import Options
from .strang_split import StrangSplit
from .wavefunction import Wavefunction


@dataclass
class _Observables:
    kinetic_energy: list[float]
    momentum: list[float]
    norm: list[float]
    population: list[float]
    position: list[float]
    potential_energy: list[float]
    total_energy: list[float]


def run(options_dict: dict):
    opt = Options(**options_dict)

    assert opt.iterations >= 0, "ITERATIONS MUST BE NON-NEGATIVE"
    assert opt.time_step > 0, "TIME STEP MUST BE POSITIVE"
    assert opt.mass > 0, "MASS MUST BE POSITIVE"
    assert opt.grid.npoint > 1, "NUMBER OF GRID POINTS MUST BE GREATER THAN 1"

    potential = opt.potential.create()
    optimized_wfns: list[Wavefunction] = []

    position_grid = generatePositionGrid(np.array(opt.grid.limits), opt.grid.npoint)
    momentum_grid = generateMomentumGrid(np.array(opt.grid.limits), opt.grid.npoint)

    npropagations = opt.imaginary.nstate if opt.imaginary else 1

    for state_idx in range(npropagations):
        wfn = Wavefunction(potential.ndim, potential.nstate, opt.grid.npoint)

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
            imaginary=opt.imaginary is not None,
        )

        with np.printoptions(formatter={"float": "{:10.4f}".format}, suppress=True):
            print(
                f"\nSTATE {state_idx} INITIAL GAMMA: "
                f"{np.array(opt.initial_conditions.gamma, dtype=float)}\n"
            )

        print(
            f"STATE {state_idx} {'IMAGINARY' if opt.imaginary else 'REAL'} TIME PROPAGATION\n"
            f"{'ITER':>5} {'KIN (Eh)':>12} {'POT (Eh)':>12} {'TOT (Eh)':>12} "
            f"{'POS (a0)':>{(11 * wfn.ndim + 1)}} {'MOM (hb/a0)':>{(11 * wfn.ndim + 1)}} "
            f"{'POPULATION':>{(11 * wfn.nstate + 1)}} {'NORM':>9} TIME"
        )

        history = {
            "kinetic_energy": [],
            "momentum": [],
            "norm": [],
            "population": [],
            "position": [],
            "potential_energy": [],
            "total_energy": [],
        }

        for i in range(opt.iterations + 1):
            start_time = time.time()

            if i: propagator.step(wfn)

            if i and opt.imaginary and state_idx > 0: wfn.projectOut(optimized_wfns)

            log_iteration = i == 0 or i == opt.iterations or (i % opt.log_interval == 0)

            if opt.write.kinetic_energy or opt.write.total_energy or log_iteration:
                kinetic_energy = wfn.kineticEnergy(momentum_grid, opt.mass)
                if opt.write.kinetic_energy:
                    history["kinetic_energy"].append(kinetic_energy)

            if opt.write.momentum or log_iteration:
                momentum = wfn.momentum(momentum_grid)
                if opt.write.momentum:
                    history["momentum"].append(momentum)

            if opt.write.norm or log_iteration:
                norm = wfn.norm()
                if opt.write.norm:
                    history["norm"].append(norm)

            if opt.write.population or log_iteration:
                population = wfn.population()
                if opt.write.population:
                    history["population"].append(population)

            if opt.write.position or log_iteration:
                position = wfn.position(position_grid)
                if opt.write.position:
                    history["position"].append(position)

            if opt.write.potential_energy or opt.write.total_energy or log_iteration:
                potential_energy = wfn.potentialEnergy(position_grid, potential)
                if opt.write.potential_energy:
                    history["potential_energy"].append(potential_energy)

            if opt.write.total_energy or log_iteration:
                total_energy = kinetic_energy + potential_energy
                if opt.write.total_energy:
                    history["total_energy"].append(total_energy)

            if not log_iteration: continue

            duration = datetime.timedelta(seconds=time.time() - start_time)

            with np.printoptions(formatter={"float": "{:10.4f}".format}, suppress=True):
                print(
                    f"{i:5d} "
                    f"{kinetic_energy:12.6f} {potential_energy:12.6f} {total_energy:12.6f} "
                    f"{position} {momentum} {population} {norm:1.3e} {duration}"
                )

        optimized_wfns.append(wfn)

        _export_history(
            history,
            opt.write,
            opt.iterations,
            opt.time_step,
            state_idx,
            npropagations > 1,
        )

def _export_history(
    history: dict[str, list],
    write_options,
    iterations: float,
    dt: float,
    state_idx: int,
    multistate: bool,
):
    for field, path in write_options:
        if path is None: continue

        times = np.arange(iterations + 1) * dt
        data = np.atleast_2d(history[field])
        output = np.column_stack((times, data))

        basename, extension = path.split(".")[0], path.split(".")[1]
        multistate_path = f"{basename}_STATE-{state_idx:02}.{extension}"
        filename = path if not multistate else multistate_path

        np.savetxt(
            filename,
            output,
            header=f"{output.shape[0]} {output.shape[1]}",
            comments="",
            fmt="%20.14f",
        )
