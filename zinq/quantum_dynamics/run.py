import datetime
import time
from dataclasses import dataclass
from typing import Any

import numpy as np

from ..potential import Potential
from .grid import generate_momentum_grid, generate_position_grid
from .options import Options, WriteOptions
from .strang_split import StrangSplit
from .wavefunction import Wavefunction


@dataclass
class _Observables:
    kinetic_energy: Any = None
    momentum: Any = None
    norm: Any = None
    population: Any = None
    position: Any = None
    potential_energy: Any = None
    total_energy: Any = None


@dataclass
class _ObservablesParams:
    wfn: Wavefunction
    potential: Potential
    position_grid: list[np.ndarray]
    momentum_grid: list[np.ndarray]
    mass: float


def run(options_dict: dict):
    opt = Options(**options_dict)

    assert opt.iterations >= 0, "ITERATIONS MUST BE NON-NEGATIVE"
    assert opt.time_step > 0, "TIME STEP MUST BE POSITIVE"
    assert opt.mass > 0, "MASS MUST BE POSITIVE"
    assert opt.grid.npoint > 1, "NUMBER OF GRID POINTS MUST BE GREATER THAN 1"

    potential = opt.potential.create()
    optimized_wfns: list[Wavefunction] = []

    position_grid = generate_position_grid(np.array(opt.grid.limits), opt.grid.npoint)
    momentum_grid = generate_momentum_grid(np.array(opt.grid.limits), opt.grid.npoint)

    npropagations = opt.imaginary.nstate if opt.imaginary else 1

    for state_idx in range(npropagations):
        wfn = Wavefunction(potential.ndim, potential.nstate, opt.grid.npoint)

        wfn.initialize_gaussian(
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

        history = _Observables(
            kinetic_energy=[],
            momentum=[],
            norm=[],
            population=[],
            position=[],
            potential_energy=[],
            total_energy=[],
        )

        obs_params = _ObservablesParams(
            wfn=wfn,
            potential=potential,
            position_grid=position_grid,
            momentum_grid=momentum_grid,
            mass=opt.mass,
        )

        for i in range(opt.iterations + 1):
            start_time = time.time()

            if i: propagator.step(wfn)

            if i and opt.imaginary and state_idx > 0: wfn.project_out(optimized_wfns)

            log_iteration = i == 0 or i == opt.iterations or (i % opt.log_interval == 0)

            obs = _calculate_observables(
                params=obs_params,
                write=opt.write,
                history=history,
                log_iteration=log_iteration,
            )

            if not log_iteration: continue

            duration = datetime.timedelta(seconds=time.time() - start_time)

            with np.printoptions(formatter={"float": "{:10.4f}".format}, suppress=True):
                print(
                    f"{i:5d} {obs.kinetic_energy:12.6f} "
                    f"{obs.potential_energy:12.6f} {obs.total_energy:12.6f} {obs.position} "
                    f"{obs.momentum} {obs.population} {obs.norm:1.3e} {duration}"
                )

        optimized_wfns.append(wfn)

        time_range = np.arange(opt.iterations + 1) * opt.time_step

        _export_history(history, opt.write, time_range, (state_idx, npropagations))


def _calculate_observables(
    params: _ObservablesParams,
    write: WriteOptions,
    history: _Observables,
    log_iteration: bool,
) -> _Observables:
    obs = _Observables()

    def _update(field, value):
        setattr(obs, field, value)
        if getattr(write, field):
            getattr(history, field).append(value)

    if write.kinetic_energy or write.total_energy or log_iteration:
        value = params.wfn.kinetic_energy(params.momentum_grid, params.mass)
        _update("kinetic_energy", value)

    if write.momentum or log_iteration:
        value = params.wfn.momentum(params.momentum_grid)
        _update("momentum", value)

    if write.norm or log_iteration:
        value = params.wfn.norm()
        _update("norm", value)

    if write.population or log_iteration:
        value = params.wfn.population()
        _update("population", value)

    if write.position or log_iteration:
        value = params.wfn.position(params.position_grid)
        _update("position", value)

    if write.potential_energy or write.total_energy or log_iteration:
        value = params.wfn.potential_energy(params.position_grid, params.potential)
        _update("potential_energy", value)

    if write.total_energy or log_iteration:
        value = obs.kinetic_energy + obs.potential_energy
        _update("total_energy", value)

    return obs


def _export_history(
    history: _Observables,
    write_options: WriteOptions,
    time_range: np.ndarray,
    istate_info: tuple[int, int],
):
    for field, path in write_options:
        if path is None: continue

        if not (data := getattr(history, field)): continue

        output = np.column_stack((time_range, np.array(data)))

        basename, extension = path.split(".")[0], path.split(".")[1]
        multistate_path = f"{basename}_STATE-{istate_info[0]:02}.{extension}"
        filename = multistate_path if istate_info[1] > 1 else path

        np.savetxt(
            filename,
            output,
            header=f"{output.shape[0]} {output.shape[1]}",
            comments="",
            fmt="%20.14f",
        )
