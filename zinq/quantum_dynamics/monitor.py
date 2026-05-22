import datetime
import time
from types import SimpleNamespace
from typing import Any

import numpy as np

from .grid import Grid
from .hamiltonian import Hamiltonian
from .options import Options
from .result import Result
from .wavefunction import Wavefunction


class Monitor:
    def __init__(self, *, idx: int, grid: Grid, H: Hamiltonian, opt: Options):
        self.idx, self.grid, self.H, self.opt = idx, grid, H, opt

        def get_history() -> dict[str, list]:
            return {k: [] for k, v in self.get_write_map(False).items() if v}

        self.history, self.latest, self.start_time = get_history(), {}, time.time()

        self._print_header()

    def export(self, dt: float) -> None:
        if self.opt.write is None: return

        for field, path in self.opt.write.model_dump().items():
            if not path: continue

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

        return {k: k in active for k in (set(Result.__annotations__) | active) - {"spectrum"}}
    
    def record(self, i: int, wfn: Wavefunction, wfn_0: Wavefunction | None, force_log: bool = False) -> None:
        wfn_adia, to_calc = wfn.to_adia(self.H) if self.opt.adiabatic else wfn, self._to_calc(i, force_log)

        ops = {
            "norm": lambda: wfn.norm(self.grid),
            "population": lambda: wfn_adia.pop(self.grid),
            "position": lambda: wfn.pos(self.grid),
            "momentum": lambda: wfn.mom(self.grid),
            "potential_energy": lambda: wfn.pe(self.grid, self.H),
            "kinetic_energy": lambda: wfn.ke(self.grid, self.H),
            "total_energy": lambda: wfn.pe(self.grid, self.H) + wfn.ke(self.grid, self.H),
            "autocorrelation": lambda: wfn_0.overlap(wfn, self.grid) if wfn_0 else None,
            "wavefunction": lambda: wfn_adia.data.copy(),
        }

        obs = {k: (ops[k]() if k in to_calc else None) for k in ops}

        self._update_history(obs)

        if self._is_log(i, force_log):
            self._print_step(SimpleNamespace(**obs), i)

    def _format_data(self, name: str, dt: float) -> np.ndarray:
        if name == "spectrum":
            return self._get_spectrum(np.array(self.history["autocorrelation"]))

        def normal_data(data) -> np.ndarray:
            return data.reshape(data.shape[0], -1)

        def normal_grid(data) -> list[np.ndarray]:
            return [np.arange(len(data)) * dt]

        def wfn_data(data) -> np.ndarray:
            return np.moveaxis(data, 0, -2).reshape(-1, data.shape[0] * data.shape[-1])

        def wfn_grid(data) -> list[np.ndarray]:
            return [r.ravel() for r in self.grid.pos]

        data, funcs = np.array(self.history[name]), {
            "wavefunction": (wfn_data, wfn_grid)
        }

        if np.iscomplexobj(data := funcs.get(name, (normal_data, normal_grid))[0](data)):
            data = np.ascontiguousarray(data).view(np.real(data).dtype)

        return np.column_stack((*funcs.get(name, (normal_data, normal_grid))[1](data), data))
    
    def _get_spectrum(self, acf: np.ndarray) -> np.ndarray:
        dt, times =self.opt.time_step, np.arange(-len(acf) + 1, len(acf)) * self.opt.time_step

        def get_omega(acf) -> np.ndarray:
            return np.fft.fftshift(2 * np.pi * np.fft.fftfreq(len(acf), dt))

        def get_spectrum(acf) -> np.ndarray:
            return np.real(np.fft.fftshift(np.fft.fft(np.fft.ifftshift(np.conj(acf))) * dt))

        def get_tau(acf) -> float:
            return -np.log(1e-4) / (len(acf) * dt)**2

        def damp_acf(acf, tau: float) -> np.ndarray:
            return acf * np.exp(-tau * times**2)

        def pad_acf(acf) -> np.ndarray:
            return np.pad(acf, 2 * [10 * len(acf)], mode="constant")

        def symmetrize_acf(acf) -> np.ndarray:
            return np.concatenate((np.conj(acf[1:][::-1]), acf))

        acf = pad_acf(damp_acf(symmetrize_acf(acf), get_tau(acf)))

        return np.column_stack((get_omega(acf), get_spectrum(acf)))

    def _is_log(self, i: int, force_log: bool = False) -> bool:
        iters, interval = self.opt.iterations, self.opt.log_interval

        return i == 0 or (iters and i == iters) or (i % interval == 0) or force_log

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

    def _to_calc(self, i: int, force_log: bool = False) -> set[str]:
        return {k for k, v in self.get_write_map(self._is_log(i, force_log)).items() if v}

    def _update_history(self, obs: dict[str, Any]) -> None:
        self.latest.update({k: v for k, v in obs.items() if v is not None})

        for k in self.history:
            if obs[k] is not None: self.history[k].append(obs[k])
