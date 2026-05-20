import datetime
import time
from itertools import count

import numpy as np

from ..potential import Potential
from .grid import Grid
from .hamiltonian import Hamiltonian
from .initial_conditions import InitialConditions
from .observables import Observables
from .options import Options
from .results import Results, State
from .strang_split import StrangSplit
from .wavefunction import Wavefunction


def run(opt: Options) -> Results:
    results, opt_wfn, nsim = [], [], opt.imaginary.nstate if opt.imaginary else 1

    for i in range(nsim):

        grid, wfn, H, pot, prop, start_time = *_init(opt), time.time()

        _print_header(i, opt.initial_conditions.gamma, grid.ndim, wfn.nstate, opt.imaginary is not None)

        obs: Observables = Observables.__new__(Observables)

        for j in range(opt.iterations + 1) if opt.iterations else count():
            if pot.is_td and j > 0:
                H._update_V(grid, pot, j * prop.dt)

            if j > 0:
                prop.step(wfn, grid, H, pot)

            if opt.imaginary and i > 0:
                wfn.project_out(opt_wfn, grid)

            log = j == 0 or j == opt.iterations or (j % opt.log_interval == 0)

            if log:
                obs = _calc_obs(wfn, grid, H, opt.adiabatic)

            if log: _, start_time = _print_step(j, obs, start_time), time.time()

        if opt.imaginary and nsim > 1 and i < nsim - 1:
            opt_wfn.append(wfn)

        results.append(State(
            total_energy=obs.e,
            kinetic_energy=obs.ke,
            potential_energy=obs.pe,
            position=obs.pos,
            momentum=obs.mom,
            norm=obs.norm,
            population=obs.pop,
        ))

    return Results(states=results)


def _calc_obs(wfn_dia: Wavefunction, grid: Grid, H: Hamiltonian, adia: bool) -> Observables:
    wfn = wfn_dia if not adia else wfn_dia.to_adia(H)

    norm = wfn.norm(grid)
    pop = wfn.pop(grid)
    pos = wfn.pos(grid)
    mom = wfn.mom(grid)

    pe = wfn_dia.pe(grid, H)
    ke = wfn_dia.ke(grid, H)

    return Observables(norm=norm, pop=pop, pos=pos, mom=mom, pe=pe, ke=ke, e=pe + ke)


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
        H=H,
        nstate=opt.hamiltonian.potential.nstate,
    )
    prop = StrangSplit(
        H=H,
        dt=opt.time_step,
        imaginary=opt.imaginary is not None,
    )
    return grid, wfn, H, opt.hamiltonian.potential, prop


def _print_header(idx: int, gamma: list[float], ndim: int, nstate: int, imag: bool):
        mode = "IMAGINARY" if imag else "REAL"

        with np.printoptions(formatter={"float": "{:10.4f}".format}, suppress=True):
            print(f"\nSTATE {idx} INITIAL GAMMA: {np.array(gamma, float)}\n")

        print(f"STATE {idx} {mode} TIME PROPAGATION")

        columns = [
            f"{'ITER':>7}",
            f"{'KIN (Eh)':>12}",
            f"{'POT (Eh)':>12}",
            f"{'TOT (Eh)':>12}",
            f"{'POS (a0)':>{11 * ndim + 1}}",
            f"{'MOM (hb/a0)':>{11 * ndim + 1}}",
            f"{'POPULATION':>{11 * nstate + 1}}",
            f"{'NORM':>9}",
            "TIME"
        ]

        print(" ".join(columns))


def _print_step(i: int, obs: Observables, time_from: float):
    with np.printoptions(formatter={"float": "{:10.4f}".format}, suppress=True):
        duration = datetime.timedelta(seconds=time.time() - time_from)

        columns = [
            f"{i:7d}",
            f"{obs.ke:12.6f}",
            f"{obs.pe:12.6f}",
            f"{obs.e:12.6f}",
            f"{obs.pos}",
            f"{obs.mom}",
            f"{obs.pop}",
            f"{obs.norm:1.3e}",
            str(duration)
        ]
        
        print(" ".join(columns), flush=True)
