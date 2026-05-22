from itertools import count
from typing import Any

import numpy as np

from .absorber import Absorber
from .grid import Grid
from .hamiltonian import Hamiltonian
from .initial_conditions import InitialConditions
from .monitor import Monitor
from .options import Options
from .result import Result
from .strang_split import StrangSplit
from .wavefunction import Wavefunction


class Runner:
    def __init__(self, opt: Options):
        def init_absorber() -> Absorber | None:
            return Absorber.from_options(opt.hamiltonian.absorber) if opt.hamiltonian.absorber else None

        def init_grid() -> Grid:
            return Grid.from_options(opt.grid)

        def init_ham(grid: Grid) -> Hamiltonian:
            return Hamiltonian(grid, opt.hamiltonian.potential, opt.hamiltonian.mass, init_absorber())

        def init_ic() -> InitialConditions:
            return InitialConditions.from_options(opt.initial_conditions)

        def init_wfn(grid: Grid, H: Hamiltonian) -> Wavefunction:
            return Wavefunction(init_ic(), grid, H)

        def init_prop(H: Hamiltonian) -> StrangSplit:
            return StrangSplit(H, opt.time_step, opt.imaginary is not None, opt.adiabatic)

        self.opt, self.grid, self.H, self.wfn, self.prop = (
            opt, (grid := init_grid()), (H := init_ham(grid)), init_wfn(grid, H), init_prop(H)
        )

    def run(self, idx: int, wfn_opt: list[Wavefunction]) -> dict[str, Any]:
        monitor = Monitor(idx=idx, grid=self.grid, H=self.H, opt=self.opt, wfn_0=self.wfn)

        decay = np.zeros(self.H.pot.nstate)

        for j in range(self.opt.iterations + 1) if self.opt.iterations else count():
            if j > 0: decay += self.prop.step(self.wfn, self.grid, self.H, (j - 0.5) * self.prop.dt)

            if self.opt.imaginary and wfn_opt:
                self.wfn.project_out(wfn_opt, self.grid)

            monitor.record(j, self.wfn, decay)

            if self.opt.iterations is None and self.wfn.norm(self.grid) < self.opt.stop_norm:
                monitor.record(j, self.wfn, decay, True); break

        monitor.export(self.prop.dt)

        return monitor.get_result_dict()


def run(opt: Options) -> Result:
    results, opt_wfn, nsim = [], [], opt.imaginary.nstate if opt.imaginary else 1

    for i in range(nsim):
        result = (runner := Runner(opt)).run(i, opt_wfn)

        if opt.imaginary and nsim > 1 and i < nsim - 1:
            opt_wfn.append(runner.wfn)

        results.append(result)

    res = Result(**{k: np.array([res[k] for res in results]) for k in Result.__annotations__})

    print()

    with np.printoptions(formatter={"float": "{:10.6f}".format}, suppress=True):
        for i, pop in enumerate(res.population): print(f"SIMULATION {i:02} POPULATIONS: {pop}")

    return res
