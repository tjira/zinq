"""Module for monitoring and exporting quantum dynamics simulation states and observables."""

import datetime
import time
from pathlib import Path
from types import SimpleNamespace
from typing import Any

import numpy as np

from .options import Options
from .result import Result
from .system import System
from .wavefunction import Wavefunction


class Monitor:
    """
    Observer and recorder for tracking simulation quantities and writing output files.

    Attributes
    ----------
    idx : int
        The simulation state index.
    system : System
        The physical system under observation.
    opt : Options
        The configuration options.
    adia : bool
        Whether variables are evaluated in the adiabatic representation.
    wfn_0 : Wavefunction | None
        Copy of the initial wavefunction for autocorrelation calculations.
    history : dict[str, list]
        Historical records of all calculated expectation values.
    latest : dict[str, Any]
        The most recent values of expectation values.
    start_time : float
        The start time of simulation.

    """

    idx: int
    system: System
    opt: Options
    adia: bool
    wfn_0: Wavefunction | None
    history: dict[str, list]
    latest: dict[str, Any]
    start_time: float

    def __init__(self, *, idx: int, system: System, opt: Options) -> None:
        """
        Initialize the Monitor.

        Parameters
        ----------
        idx : int
            The state index.
        system : System
            The physical system.
        opt : Options
            The simulation options.

        """
        self.idx, self.system, self.opt, self.adia = idx, system, opt, opt.adiabatic

        self.wfn_0 = system.wfn.copy() if "autocorrelation" in self.get_write_map() else None

        self.history = {k: [] for k in self.get_write_map(log=False)}
        self.latest = {}
        self.start_time = time.time()

        self._print_header()

    def export(self, dt: float) -> None:
        """
        Export recorded observables to output files as configured in write options.

        Parameters
        ----------
        dt : float
            The simulation time step.

        """
        if self.opt.write is None:
            return

        nsim = self.opt.imaginary.nstate if self.opt.imaginary else 1

        for field, path in self.opt.write.model_dump().items():
            if not path:
                continue

            output_path = path
            if nsim > 1:
                p = Path(path)
                output_path = str(p.parent / f"{p.stem}_STATE-{self.idx:02}{p.suffix}")

            data = self._format_data(field, dt)
            header = f"{data.shape[0]} {data.shape[1]}"

            np.savetxt(output_path, data, header=header, comments="", fmt="%20.14f")

    def get_result_dict(self) -> dict[str, Any]:
        """
        Get a dictionary containing the latest recorded values.

        Returns
        -------
        dict[str, Any]
            The latest recorded observables.

        """
        return {k: v for k, v in self.latest.items() if v is not None}

    def get_write_map(self, *, log: bool = False) -> set[str]:
        """
        Get the set of observables that need to be calculated and recorded.

        Parameters
        ----------
        log : bool
            Whether to include observables required for console logging.

        Returns
        -------
        set[str]
            Set of active observable names.

        """
        active = set(self.opt.write.model_dump(exclude_none=True)) if self.opt.write else set()

        if "spectrum" in active:
            active.add("autocorrelation")

        if log:
            active |= {
                "potential_energy",
                "kinetic_energy",
                "total_energy",
                "position",
                "momentum",
                "population",
                "norm",
            }

        return (active | (set(Result.__annotations__) & active)) - {"spectrum"}

    def record(self, i: int, decay: np.ndarray, *, force_log: bool = False) -> None:
        """
        Record the state of the system at time step i.

        Parameters
        ----------
        i : int
            The current time step.
        decay : np.ndarray
            The accumulated norm decay.
        force_log : bool
            Whether to force logging of this step.

        """
        if not (to_calc := self._to_calc(i, force_log=force_log)):
            return

        wfn_adia = self.system.wfn.to_adia(self.system.ham) if self.adia else self.system.wfn

        ops = {
            "norm": lambda: self.system.norm,
            "population": lambda: wfn_adia.pop(self.system.grid),
            "position": lambda: self.system.pos,
            "momentum": lambda: self.system.mom,
            "potential_energy": lambda: self.system.pe,
            "kinetic_energy": lambda: self.system.ke,
            "total_energy": lambda: self.system.pe + self.system.ke,
            "autocorrelation": lambda: (
                self.wfn_0.overlap(self.system.wfn, self.system.grid)
                if self.wfn_0 is not None
                else None
            ),
            "wavefunction": wfn_adia.data.copy,
        }

        obs = {k: (ops[k]() if k in to_calc else None) for k in ops}

        if self.system.ham.absorber and obs["population"] is not None:
            obs["population"] += decay

        self._update_history(obs)

        if self._is_log(i, force_log=force_log):
            self._print_step(SimpleNamespace(**obs), i)

    def _normal_data(self, data: np.ndarray) -> np.ndarray:
        """Reshape history data to 2D shape."""
        return data.reshape(data.shape[0], -1)

    def _normal_grid(self, data: np.ndarray, dt: float) -> list[np.ndarray]:
        """Generate positional/temporal grid coordinates for normal data."""
        return [np.arange(len(data)) * dt]

    def _wfn_data(self, data: np.ndarray) -> np.ndarray:
        """Rearrange complexity and spatial dims of wavefunction history."""
        return np.moveaxis(data, 0, -2).reshape(-1, data.shape[0] * data.shape[-1])

    def _wfn_grid(self, _data: np.ndarray) -> list[np.ndarray]:
        """Retrieve flattened grid positions."""
        return [r.ravel() for r in self.system.grid.pos]

    def _format_data(self, name: str, dt: float) -> np.ndarray:
        """Format the history of an observable into a plottable/exportable array."""
        if name == "spectrum":
            return self._get_spectrum(np.array(self.history["autocorrelation"]))

        funcs = {
            "wavefunction": (
                self._wfn_data,
                self._wfn_grid,
            )
        }
        data = np.array(self.history[name])

        get_data_func, get_grid_func = funcs.get(
            name,
            (
                self._normal_data,
                lambda d: self._normal_grid(d, dt),
            ),
        )

        data = get_data_func(data)

        if np.iscomplexobj(data):
            data = np.ascontiguousarray(data).view(np.real(data).dtype)

        return np.column_stack((*get_grid_func(data), data))

    def _get_spectrum(self, acf: np.ndarray) -> np.ndarray:
        """Compute the power spectrum from the autocorrelation function."""
        dt = self.opt.time_step
        n = len(acf)
        times = np.arange(-n + 1, n) * dt

        tau = -np.log(1e-4) / (n * dt) ** 2
        acf = np.concatenate((np.conj(acf[1:][::-1]), acf))
        acf = acf * np.exp(-tau * times**2)

        pad_width = 10 * len(acf)
        acf = np.pad(acf, (pad_width, pad_width), mode="constant")

        omega = np.fft.fftshift(2 * np.pi * np.fft.fftfreq(len(acf), dt))
        spectrum = np.real(np.fft.fftshift(np.fft.fft(np.fft.ifftshift(np.conj(acf))) * dt))

        return np.column_stack((omega, spectrum))

    def _is_log(self, i: int, *, force_log: bool = False) -> bool:
        """Determine if logging should occur for time step i."""
        iters, interval = self.opt.iterations, self.opt.log_interval

        return i == 0 or (iters and i == iters) or (i % interval == 0) or force_log

    def _print_header(self) -> None:
        """Print the console logging table header."""
        mode = "IMAGINARY" if self.opt.imaginary else "REAL"
        g = self.opt.initial_conditions.gamma

        with np.printoptions(formatter={"float": "{:10.4f}".format}, suppress=True):
            print(f"\nSTATE {self.idx} INITIAL GAMMA: {np.array(g, float)}\n")

        print(f"STATE {self.idx} {mode} TIME PROPAGATION")

        cols = [
            f"{'ITER':>7}",
            f"{'KIN (Eh)':>12}",
            f"{'POT (Eh)':>12}",
            f"{'TOT (Eh)':>12}",
            f"{'POS (a0)':>{11 * self.system.grid.ndim + 1}}",
            f"{'MOM (hb/a0)':>{11 * self.system.grid.ndim + 1}}",
            f"{'POPULATION':>{11 * self.system.ham.pot.nstate + 1}}",
            f"{'NORM':>9}",
            "TIME",
        ]

        print(" ".join(cols))

    def _print_step(self, obs: SimpleNamespace, i: int) -> None:
        """Print physical quantities for a single time step."""
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
                str(duration),
            ]

            print(" ".join(cols), flush=True)

            self.start_time = time.time()

    def _to_calc(self, i: int, *, force_log: bool = False) -> set[str]:
        """Get the observables to calculate for this step."""
        return self.get_write_map(log=self._is_log(i, force_log=force_log))

    def _update_history(self, obs: dict[str, Any]) -> None:
        """Update historical tracking and latest state of observables."""
        self.latest.update({k: v for k, v in obs.items() if v is not None})

        for k in self.history:
            if obs[k] is not None:
                self.history[k].append(obs[k])
