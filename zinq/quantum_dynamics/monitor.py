import datetime
import time
from types import SimpleNamespace
from typing import Any

import numpy as np

from .grid import Grid
from .hamiltonian import Hamiltonian
from .options import Options, WriteConfig
from .result import Result
from .wavefunction import Wavefunction


class Monitor:
    def __init__(self, *, idx: int, grid: Grid, H: Hamiltonian, opt: Options):
        self.__dict__.update({k: v for k, v in locals().items() if k != "self"})

        def get_history() -> dict[str, list[Any]]:
            return {k: [] for k, v in self.get_write_map(True).items() if v}

        self.history, self.start_time, _ = get_history(), time.time(), self._print_header()

    @property
    def latest(self) -> Any:
        return {k: v[-1] if k in self.history and v else None for k, v in self.history.items()}

    def record(self, j: int, wfn: Wavefunction, wfn_0: Wavefunction | None) -> None:
        log = j == 0 or (self.opt.iterations and j == self.opt.iterations) or (j % self.opt.log_interval == 0)

        obs_to_calc = {k for k, v in self.get_write_map(log).items() if v}

        wfn_adia = wfn.to_adia(self.H) if self.opt.adiabatic else wfn

        pe = wfn.pe(self.grid, self.H) if "potential_energy" in obs_to_calc or "total_energy" in obs_to_calc else None
        ke = wfn.ke(self.grid, self.H) if "kinetic_energy" in obs_to_calc or "total_energy" in obs_to_calc else None

        obs = {
            "norm": wfn.norm(self.grid) if "norm" in obs_to_calc else None,
            "population": wfn_adia.pop(self.grid) if "population" in obs_to_calc else None,
            "position": wfn.pos(self.grid) if "position" in obs_to_calc else None,
            "momentum": wfn.mom(self.grid) if "momentum" in obs_to_calc else None,
            "potential_energy": pe,
            "kinetic_energy": ke,
            "total_energy": pe + ke if "total_energy" in obs_to_calc and pe is not None and ke is not None else None,
            "autocorrelation": wfn_0.overlap(wfn, self.grid) if "autocorrelation" in obs_to_calc and wfn_0 else None,
            "wavefunction": wfn_adia.data.copy() if "wavefunction" in obs_to_calc else None,
        }

        for k, v in (kv for kv in obs.items() if kv[1] is not None): self.history[k].append(v)

        if log: self._print_step(SimpleNamespace(**obs), j)

    def export(self, dt: float) -> None:
        if self.opt.write is None: return

        for field, path in self.opt.write.model_dump().items():
            if not path or field == "spectrum": continue

            data = self._format_data(field, dt); header = f"{data.shape[0]} {data.shape[1]}"

            np.savetxt(path, data, header=header, comments="", fmt="%20.14f")

    def get_result_dict(self) -> dict[str, Any]:
        return {k: v for k, v in self.latest.items() if v is not None}

    def get_write_map(self, log: bool = False) -> dict[str, bool]:
        active = set(self.opt.write.model_dump(exclude_none=True)) if self.opt.write else set()

        if "spectrum" in active: active.add("autocorrelation")

        if log: active |= {
            "potential_energy",
            "kinetic_energy",
            "total_energy",
            "position",
            "momentum",
            "population",
            "norm"
        }

        return {k: k in active for k in set(Result.__annotations__) | active}

    def _format_data(self, name: str, dt: float) -> np.ndarray:
        def normal_data(data) -> np.ndarray:
            return data.reshape(data.shape[0], -1)

        def normal_grid(data) -> np.ndarray:
            return [np.arange(len(data)) * dt]

        def wfn_data(data) -> np.ndarray:
            return np.moveaxis(data, 0, -2).reshape(-1, data.shape[0] * data.shape[-1])

        def wfn_grid(data) -> np.ndarray:
            return [r.ravel() for r in self.grid.pos]

        data, funcs = np.array(self.history[name]), {
            "wavefunction": (wfn_data, wfn_grid)
        }

        if np.iscomplexobj(data := funcs.get(name, (normal_data, normal_grid))[0](data)):
            data = np.ascontiguousarray(data).view(np.real(data).dtype)

        return np.column_stack((*funcs.get(name, (normal_data, normal_grid))[1](data), data))

    def _print_header(self) -> None:
        mode, g = "IMAGINARY" if self.opt.imaginary else "REAL", self.opt.initial_conditions.gamma

        with np.printoptions(formatter={"float": "{:10.4f}".format}, suppress=True):
            print(f"\nSTATE {self.idx} INITIAL GAMMA: {np.array(g, float)}\n")

        print(f"STATE {self.idx} {mode} TIME PROPAGATION")

        cols = [
            f"{'ITER':>7}", f"{'KIN (Eh)':>12}", f"{'POT (Eh)':>12}", f"{'TOT (Eh)':>12}",
            f"{'POS (a0)':>{11 * self.grid.ndim + 1}}",
            f"{'MOM (hb/a0)':>{11 * self.grid.ndim + 1}}",
            f"{'POPULATION':>{11 * self.H.pot.nstate + 1}}",
            f"{'NORM':>9}", "TIME"
        ]

        print(" ".join(cols))

    def _print_step(self, obs: SimpleNamespace, i: int) -> None:
        with np.printoptions(formatter={"float": "{:10.4f}".format}, suppress=True):
            duration = datetime.timedelta(seconds=time.time() - self.start_time)

            cols = [
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

            print(" ".join(cols), flush=True)

            self.start_time = time.time()
