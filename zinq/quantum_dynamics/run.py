import datetime
import time
from dataclasses import dataclass
from itertools import count
from types import SimpleNamespace
from typing import Any

import numpy as np

from .grid import Grid
from .hamiltonian import Hamiltonian
from .history import History
from .initial_conditions import InitialConditions
from .options import HamiltonianConfig, Options, WriteConfig
from .result import Result
from .strang_split import StrangSplit
from .wavefunction import Wavefunction


@dataclass(kw_only=True)
class Runner:
    H: Hamiltonian
    grid: Grid
    opt: Options
    prop: StrangSplit
    wfn: Wavefunction

    def run(self, idx: int, wfn_opt: list[Wavefunction]) -> dict[str, Any]:
        history, start_time = History(), time.time()

        _print_header(idx, self.opt.initial_conditions.gamma, self.wfn, self.prop.imag)

        wfn_0 = self.wfn.copy() if _obs_map(self.opt.write)["autocorrelation"] else None

        for j in range(self.opt.iterations + 1) if self.opt.iterations else count():
            if j > 0: self.prop.step(self.wfn, self.grid, self.H, j * self.prop.dt)

            if self.opt.imaginary and wfn_opt:
                self.wfn.project_out(wfn_opt, self.grid)

            log = j == 0 or j == self.opt.iterations or (j % self.opt.log_interval == 0)

            history.record(**self._calc_obs(_obs_map(self.opt.write, log), wfn_0))

            if log:
                _, start_time = _print_step(j, history.latest, start_time), time.time()
        
        self._export_obs(history, self.opt.write, self.prop.dt)

        return {obs: getattr(history.latest, obs) for obs in Result.__annotations__}

    def _calc_obs(self, write_map: dict[str, bool], wfn_0: Wavefunction | None) -> dict[str, Any]:
        wfn = self.wfn.to_adia(self.H) if self.opt.adiabatic else self.wfn

        pe = self.wfn.pe(self.grid, self.H) if write_map["potential_energy"] or write_map["total_energy"] else None
        ke = self.wfn.ke(self.grid, self.H) if write_map["kinetic_energy"] or write_map["total_energy"] else None

        return {
            "norm": self.wfn.norm(self.grid) if write_map["norm"] else None,
            "population": wfn.pop(self.grid) if write_map["population"] else None,
            "position": self.wfn.pos(self.grid) if write_map["position"] else None,
            "momentum": self.wfn.mom(self.grid) if write_map["momentum"] else None,
            "potential_energy": pe,
            "kinetic_energy": ke,
            "total_energy": pe + ke if write_map["total_energy"] and pe is not None and ke is not None else None,
            "autocorrelation": wfn_0.overlap(self.wfn, self.grid) if write_map["autocorrelation"] and wfn_0 else None,
            "wavefunction": wfn.data.copy() if write_map["wavefunction"] else None,
        }


    def _export_obs(self, history: History, write_opts: WriteConfig | None, dt: float) -> None:
        if write_opts is None: return

        for field, path in write_opts.model_dump().items():
            if not path or field == "spectrum": continue

            data = history.get(field, self.grid, dt)

            np.savetxt(path, data, header=f"{data.shape[0]} {data.shape[1]}", comments="", fmt="%20.14f")


def run(opt: Options) -> Result:
    states, opt_wfn, nsim = [], [], opt.imaginary.nstate if opt.imaginary else 1

    for i in range(nsim):
        state = (runner := Runner(**_init(opt), opt=opt)).run(i, opt_wfn)

        if opt.imaginary and nsim > 1 and i < nsim - 1:
            opt_wfn.append(runner.wfn)

        states.append(state)

    return Result(**{obs: np.array([state[obs] for state in states]) for obs in Result.__annotations__})


def _init(opt: Options) -> dict:
    def _init_grid(opt: Options) -> Grid:
        return Grid.from_options(opt.grid)

    def _init_ham(opt: HamiltonianConfig, grid: Grid) -> Hamiltonian:
        return Hamiltonian(grid=grid, pot=opt.potential, m=opt.mass)

    def _init_ic(opt: Options) -> InitialConditions:
        return InitialConditions.from_options(opt.initial_conditions)

    def _init_prop(opt: Options, H: Hamiltonian) -> StrangSplit:
        return StrangSplit(H=H, dt=opt.time_step, imag=opt.imaginary is not None)

    def _init_wfn(opt: Options, grid: Grid, H: Hamiltonian) -> Wavefunction:
        return Wavefunction(ic=_init_ic(opt), grid=grid, H=H)

    grid, H = (grid := _init_grid(opt)), _init_ham(opt.hamiltonian, grid)

    return {"grid": grid, "H": H, "wfn": _init_wfn(opt, grid, H), "prop": _init_prop(opt, H)}


def _obs_map(write_opts: WriteConfig | None, log: bool = False) -> dict[str, bool]:
    def is_set(attr: str) -> bool:
        return write_opts is not None and getattr(write_opts, attr) is not None

    spectrum = is_set("spectrum")

    return {
        "autocorrelation": is_set("autocorrelation") or spectrum,
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


def _print_header(idx: int, gamma: list[float], wfn: Wavefunction, imag: bool) -> None:
    mode = "IMAGINARY" if imag else "REAL"

    with np.printoptions(formatter={"float": "{:10.4f}".format}, suppress=True):
        print(f"\nSTATE {idx} INITIAL GAMMA: {np.array(gamma, float)}\n")

    print(f"STATE {idx} {mode} TIME PROPAGATION")

    columns = [
        f"{'ITER':>7}",
        f"{'KIN (Eh)':>12}",
        f"{'POT (Eh)':>12}",
        f"{'TOT (Eh)':>12}",
        f"{'POS (a0)':>{11 * wfn.ndim + 1}}",
        f"{'MOM (hb/a0)':>{11 * wfn.ndim + 1}}",
        f"{'POPULATION':>{11 * wfn.nstate + 1}}",
        f"{'NORM':>9}",
        "TIME"
    ]

    print(" ".join(columns))


def _print_step(i: int, obs: SimpleNamespace, time_from: float) -> None:
    with np.printoptions(formatter={"float": "{:10.4f}".format}, suppress=True):
        duration = datetime.timedelta(seconds=time.time() - time_from)

        columns = [
            f"{i:7d}",
            f"{obs.kinetic_energy:12.6f}",
            f"{obs.potential_energy:12.6f}",
            f"{obs.total_energy:12.6f}",
            f"{obs.position}",
            f"{obs.momentum}",
            f"{obs.population}",
            f"{obs.norm:1.3e}",
            str(duration)
        ]

        print(" ".join(columns), flush=True)
