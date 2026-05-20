import datetime
import time
from itertools import count
from typing import cast

import numpy as np

from ..potential import Potential
from .grid import Grid
from .hamiltonian import Hamiltonian
from .initial_conditions import InitialConditions
from .observables import Observables
from .options import Options, WriteOptions
from .results import Results, State
from .strang_split import StrangSplit
from .wavefunction import Wavefunction


class Runner:
    H: Hamiltonian
    grid: Grid
    opt: Options
    pot: Potential
    prop: StrangSplit
    wfn: Wavefunction

    def __init__(self, opt: Options) -> None:
        self.opt, self.pot = opt, opt.hamiltonian.potential

        self.grid = Grid(
            limits=np.array(opt.grid.limits),
            npoint=opt.grid.npoint
        )
        self.H = Hamiltonian(
            grid=self.grid,
            pot=self.pot,
            m=opt.hamiltonian.mass
        )
        self.wfn = Wavefunction(
            ic=InitialConditions(
                pos=np.array(opt.initial_conditions.position),
                mom=np.array(opt.initial_conditions.momentum),
                gamma=np.array(opt.initial_conditions.gamma),
                state=opt.initial_conditions.state,
                adia=opt.initial_conditions.adiabatic,
            ),
            grid=self.grid,
            H=self.H,
            nstate=self.pot.nstate
        )
        self.prop = StrangSplit(
            H=self.H,
            dt=opt.time_step,
            imag=opt.imaginary is not None
        )

    def run(self, idx: int = 0, prev_states: list[Wavefunction] | None = None) -> State:
        obs, start_time = Observables.__new__(Observables), time.time()

        _print_header(idx, self.opt.initial_conditions.gamma, self.grid.ndim, self.wfn.nstate, self.opt.imaginary is not None)

        wfn_0 = self.wfn.copy() if _obs_map(self.opt.write)["autocorrelation"] else None

        for j in range(self.opt.iterations + 1) if self.opt.iterations else count():
            if self.pot.is_td and j > 0:
                self.H._update_V(self.grid, self.pot, j * self.prop.dt)

            if j > 0:
                self.prop.step(self.wfn, self.grid, self.H, self.pot)

            if self.opt.imaginary and prev_states:
                self.wfn.project_out(prev_states, self.grid)

            log = j == 0 or j == self.opt.iterations or (j % self.opt.log_interval == 0)

            obs = self._calc_obs(_obs_map(self.opt.write, log), wfn_0)

            if log:
                _, start_time = _print_step(j, obs, start_time), time.time()

        return State(
            total_energy=cast(float, obs.e),
            kinetic_energy=cast(float, obs.ke),
            potential_energy=cast(float, obs.pe),
            position=cast(np.ndarray, obs.pos),
            momentum=cast(np.ndarray, obs.mom),
            norm=cast(float, obs.norm),
            population=cast(np.ndarray, obs.pop),
        )

    def _calc_obs(self, write_map: dict[str, bool], wfn_0: Wavefunction | None) -> Observables:
        wfn = self.wfn.to_adia(self.H) if self.opt.adiabatic else self.wfn

        pe = self.wfn.pe(self.grid, self.H) if write_map.get("potential_energy") or write_map.get("total_energy") else None
        ke = self.wfn.ke(self.grid, self.H) if write_map.get("kinetic_energy") or write_map.get("total_energy") else None

        return Observables(
            norm=wfn.norm(self.grid) if write_map.get("norm") else None,
            pop=wfn.pop(self.grid) if write_map.get("population") else None,
            pos=wfn.pos(self.grid) if write_map.get("position") else None,
            mom=wfn.mom(self.grid) if write_map.get("momentum") else None,
            pe=pe, ke=ke, e=pe + ke if write_map.get("total_energy") and pe is not None and ke is not None else None,
            acf=wfn_0.overlap(wfn, self.grid) if write_map.get("autocorrelation") and wfn_0 else None,
            wfn=wfn.data.copy() if write_map.get("wavefunction") else None
        )


def run(opt: Options) -> Results:
    results, opt_wfn, nsim = [], [], opt.imaginary.nstate if opt.imaginary else 1

    for i in range(nsim):
        state = (runner := Runner(opt)).run(i, opt_wfn)

        if opt.imaginary and nsim > 1 and i < nsim - 1:
            opt_wfn.append(runner.wfn)

        results.append(state)

    return Results(states=results)


def _obs_map(write_opts: WriteOptions | None, log: bool = False) -> dict[str, bool]:
    def is_set(attr: str) -> bool:
        return write_opts is not None and getattr(write_opts, attr) is not None

    spectrum = is_set("spectrum")

    return {
        "autocorrelation": is_set("autocorrelation") or spectrum,
        "final_wavefunction": is_set("final_wavefunction"),
        "kinetic_energy": is_set("kinetic_energy") or log,
        "momentum": is_set("momentum") or log,
        "norm": is_set("norm") or log,
        "population": is_set("population") or log,
        "position": is_set("position") or log,
        "potential_energy": is_set("potential_energy") or log,
        "spectrum": spectrum,
        "total_energy": is_set("total_energy") or log,
        "wavefunction": is_set("wavefunction"),
    }


def _print_header(idx: int, gamma: list[float], ndim: int, nstate: int, imag: bool) -> None:
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


def _print_step(i: int, obs: Observables, time_from: float) -> None:
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
