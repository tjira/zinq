"""Module executing quantum dynamics propagation simulations."""

from itertools import count
from typing import Any

import numpy as np

from .monitor import Monitor
from .options import Options
from .result import Result
from .strang_split import StrangSplit
from .system import System
from .wavefunction import Wavefunction


class Runner:
    """
    Simulation runner for managing quantum dynamics propagation.

    Attributes
    ----------
    opt : Options
        The configuration options.
    system : System
        The physical system undergoing propagation.
    prop : StrangSplit
        The Strang split-operator propagator.

    """

    opt: Options
    system: System
    prop: StrangSplit

    def __init__(self, opt: Options) -> None:
        """
        Initialize the runner.

        Parameters
        ----------
        opt : Options
            The simulation options.

        """
        self.opt, self.system = opt, System.from_options(opt)

        self.prop = StrangSplit(
            ham=self.system.ham,
            dt=opt.time_step,
            imag=opt.imaginary is not None,
            adia=opt.adiabatic,
        )

    def run(self, idx: int, wfn_opt: list[Wavefunction]) -> dict[str, Any]:
        """
        Run the propagation simulation.

        Parameters
        ----------
        idx : int
            The simulation index (e.g. state index in imaginary time).
        wfn_opt : list[Wavefunction]
            Optimized lower-energy wavefunctions to project out.

        Returns
        -------
        dict[str, Any]
            Dictionary containing recorded expectation values over time.

        """
        monitor = Monitor(idx=idx, system=self.system, opt=self.opt)

        decay = np.zeros(self.system.ham.pot.nstate)

        for j in range(self.opt.iterations + 1) if self.opt.iterations else count():
            if j > 0:
                decay += self.prop.step(self.system, (j - 0.5) * self.prop.dt)

            if self.opt.imaginary and wfn_opt:
                self.system.wfn.project_out(wfn_opt, self.system.grid)

            monitor.record(j, decay)

            if self.opt.iterations is None and self.system.norm < self.opt.stop_norm:
                break

        monitor.export(self.prop.dt)

        return monitor.get_result_dict()


def run(opt: Options) -> Result:
    """
    Run one or multiple quantum dynamics simulations as configured by options.

    Parameters
    ----------
    opt : Options
        The configuration options.

    Returns
    -------
    Result
        The consolidated results of the simulation.

    """
    results, opt_wfn, nsim = [], [], opt.imaginary.nstate if opt.imaginary else 1

    for i in range(nsim):
        result = (runner := Runner(opt)).run(i, opt_wfn)

        if opt.imaginary and nsim > 1 and i < nsim - 1:
            opt_wfn.append(runner.system.wfn)

        results.append(result)

    return Result(**{k: np.array([res[k] for res in results]) for k in Result.__annotations__})
