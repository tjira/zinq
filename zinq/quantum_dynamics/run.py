import datetime
import time

import numpy as np

from ..potential import Potential
from .grid import generate_momentum_grid, generate_position_grid
from .options import Options
from .strang_split import StrangSplit
from .wavefunction import Wavefunction


def run(options_dict: dict):
    opt = Options(**options_dict)
    _validate(opt)
    Runner(opt).run()


def _validate(opt: Options):
    assert opt.iterations >= 0, "ITERATIONS MUST BE NON-NEGATIVE"
    assert opt.time_step > 0, "TIME STEP MUST BE POSITIVE"
    assert opt.mass > 0, "MASS MUST BE POSITIVE"
    assert opt.grid.npoint > 1, "NUMBER OF GRID POINTS MUST BE GREATER THAN 1"


class Runner:
    opt: Options
    pot: Potential
    pos_grid: list[np.ndarray]
    mom_grid: list[np.ndarray]
    optimized: list[Wavefunction] = []

    def __init__(self, opt: Options):
        self.opt, self.pot = opt, opt.potential.create()

        self.pos_grid = generate_position_grid(np.array(opt.grid.limits), opt.grid.npoint)
        self.mom_grid = generate_momentum_grid(np.array(opt.grid.limits), opt.grid.npoint)

    def run(self):
        n = self.opt.imaginary.nstate if self.opt.imaginary else 1

        for i in range(n):
            wfn = self._init_wfn()
            self._propagate(wfn, i, n)
            self.optimized.append(wfn)

    def _init_wfn(self) -> Wavefunction:
        wfn = Wavefunction(self.pot.ndim, self.pot.nstate, self.opt.grid.npoint)

        wfn.initialize_gaussian(
            self.pos_grid,
            np.array(self.opt.initial_conditions.position),
            np.array(self.opt.initial_conditions.momentum),
            np.array(self.opt.initial_conditions.gamma),
            self.opt.initial_conditions.state,
        )

        return wfn

    def _propagate(self, wfn: Wavefunction, idx: int, nstate: int):
        prop = StrangSplit(
            self.pos_grid,
            self.mom_grid,
            self.pot,
            self.opt.mass,
            self.opt.time_step,
            self.opt.imaginary is not None
        )

        history = {f: [] for f, p in self.opt.write if p is not None}

        self._print_header(idx, wfn)

        for i in range(self.opt.iterations + 1):
            start = time.time()

            if i: prop.step(wfn)

            if i and self.opt.imaginary and idx > 0: wfn.project_out(self.optimized)

            is_log = i == 0 or i == self.opt.iterations or (i % self.opt.log_interval == 0)

            obs = self._calc_obs(wfn, is_log)

            for f, v in obs.items():
                if f in history: history[f].append(v)

            if is_log:
                self._print_step(i, obs, datetime.timedelta(seconds=time.time() - start))

        self._save(history, idx, nstate)

    def _calc_obs(self, wfn: Wavefunction, is_log: bool) -> dict:
        results = {}

        def should(f): return is_log or getattr(self.opt.write, f)

        if should("kinetic_energy") or should("total_energy"):
            results["kinetic_energy"] = wfn.kinetic_energy(self.mom_grid, self.opt.mass)
        if should("momentum"):
            results["momentum"] = wfn.momentum(self.mom_grid)
        if should("norm"):
            results["norm"] = wfn.norm()
        if should("population"):
            results["population"] = wfn.population()
        if should("position"):
            results["position"] = wfn.position(self.pos_grid)
        if should("potential_energy") or should("total_energy"):
            results["potential_energy"] = wfn.potential_energy(self.pos_grid, self.pot)
        if should("total_energy"):
            results["total_energy"] = results["kinetic_energy"] + results["potential_energy"]

        return results

    def _print_header(self, idx: int, wfn: Wavefunction):
        ic, mode = self.opt.initial_conditions, "IMAGINARY" if self.opt.imaginary else "REAL"

        with np.printoptions(formatter={"float": "{:10.4f}".format}, suppress=True):
            print(f"\nSTATE {idx} INITIAL GAMMA: {np.array(ic.gamma, float)}\n")

        print(f"STATE {idx} {mode} TIME PROPAGATION")

        p_w, s_w = 11 * wfn.ndim + 1, 11 * wfn.nstate + 1

        print(
            f"{'ITER':>5} {'KIN (Eh)':>12} {'POT (Eh)':>12} {'TOT (Eh)':>12} "
            f"{'POS (a0)':>{p_w}} {'MOM (hb/a0)':>{p_w}} "
            f"{'POPULATION':>{s_w}} {'NORM':>9} TIME"
        )

    def _print_step(self, i: int, obs: dict, duration: datetime.timedelta):
        with np.printoptions(formatter={"float": "{:10.4f}".format}, suppress=True):
            print(
                f"{i:5d} {obs['kinetic_energy']:12.6f} "
                f"{obs['potential_energy']:12.6f} {obs['total_energy']:12.6f} "
                f"{obs['position']} {obs['momentum']} {obs['population']} "
                f"{obs['norm']:1.3e} {duration}"
            )

    def _save(self, history: dict, idx: int, nstate: int):
        times = np.arange(self.opt.iterations + 1) * self.opt.time_step

        for field, path in self.opt.write:

            if not path or field not in history: continue

            base, ext = path.rsplit(".", 1)

            data = np.column_stack((times, np.array(history[field])))
            fname = path if nstate == 1 else f"{base}_STATE-{idx:02}.{ext}"

            np.savetxt(
                fname,
                data,
                header=f"{data.shape[0]} {data.shape[1]}",
                comments="",
                fmt="%20.14f"
            )
