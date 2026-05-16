import datetime
import time
from itertools import count
from typing import Optional

from ..backend import np
from .grid import Grid
from .hamiltonian import Hamiltonian
from .options import Options
from .results import RunResult, StateResult
from .strang_split import StrangSplit
from .wavefunction import Wavefunction


def run(options_dict: dict) -> RunResult:
    opt = Options(**options_dict)
    _validate(opt)
    return Runner(opt).run()


def _validate(opt: Options):
    if opt.iterations: assert opt.iterations >= 0, "ITERATIONS MUST BE NON-NEGATIVE"
    assert opt.time_step > 0, "TIME STEP MUST BE POSITIVE"
    assert opt.mass > 0, "MASS MUST BE POSITIVE"
    assert opt.grid.npoint > 1, "NUMBER OF GRID POINTS MUST BE GREATER THAN 1"


class Runner:
    opt: Options
    ham: Hamiltonian
    grid: Grid
    optimized: list[Wavefunction] = []

    def __init__(self, opt: Options):
        self.opt = opt
        self.ham = Hamiltonian(opt.potential.create(), opt.mass, opt.absorbing_potential)
        self.grid = Grid(np.array(opt.grid.limits), opt.grid.npoint)

    def run(self) -> RunResult:
        n, results = self.opt.imaginary.nstate if self.opt.imaginary else 1, []

        for i in range(n):
            wfn = self._init()
            results.append(self._prop(wfn, i, n))
            self.optimized.append(wfn)

        return RunResult(states=results)

    def _init(self) -> Wavefunction:
        ndim, nstate = self.ham.potential.ndim, self.ham.potential.nstate

        pos = np.array(self.opt.initial_conditions.position)
        mom = np.array(self.opt.initial_conditions.momentum)
        gamma = np.array(self.opt.initial_conditions.gamma)
        state = self.opt.initial_conditions.state

        wfn = Wavefunction(ndim, nstate, self.opt.grid.npoint)

        wfn.init_gauss(self.grid, pos, mom, gamma, state)

        if self.opt.initial_conditions.adiabatic:
            wfn = wfn.to_diabatic(self.grid, self.ham)

        return wfn

    def _prop(self, wfn: Wavefunction, idx: int, nstate: int) -> StateResult:
        img, dt, pop_decay = self.opt.imaginary is not None, self.opt.time_step, np.zeros(wfn.nstate)

        history, iters = {f: [] for f, p in self.opt.write if p is not None}, 0

        need_acf = bool(self.opt.write.autocorrelation) or bool(self.opt.write.spectrum)

        if need_acf and "autocorrelation" not in history: history["autocorrelation"] = []

        self._head(idx, wfn)

        propagator = StrangSplit(self.grid, self.ham, dt, img, self.opt.adiabatic)

        wfn_0 = Wavefunction.from_data(wfn.data.copy(), wfn.measure) if need_acf else None

        for i in range(self.opt.iterations + 1) if self.opt.iterations else count():
            start, current_time, iters = time.time(), i * dt, iters + 1

            if i: pop_decay += propagator.step(wfn, self.grid, current_time - dt)

            if i and self.opt.imaginary and idx > 0: wfn.project_out(self.optimized)

            log = i == 0 or i == self.opt.iterations or (i % self.opt.log_interval == 0)

            obs = self._obs(wfn, log, pop_decay, current_time, wfn_0)

            for f, v in obs.items():
                if f in history: history[f].append(v)

            elapsed = datetime.timedelta(seconds=time.time() - start)

            if log: self._log_step(i, obs, elapsed)

            if self.opt.absorbing_potential and obs["norm"] < self.opt.absorbing_potential.stop_norm:
                if not log: self._log_step(i, self._obs(wfn, True, pop_decay, current_time), elapsed)
                print(f"\nCAP STOP NORM REACHED, STOPPING PROPAGATION")
                break

        if "final_wavefunction" in history:
            wfn_final = wfn.to_adiabatic(self.grid, self.ham, iters * dt) if self.opt.adiabatic else wfn
            history["final_wavefunction"] = [wfn_final.data.copy()]

        if "spectrum" in history:
            history["spectrum"] = [self._get_spectrum(np.array(history["autocorrelation"]))]

        final_obs = self._obs(wfn, True, pop_decay, iters * dt)

        with np.printoptions(formatter={"float": "{:10.4f}".format}, suppress=True):
            print(f"\nFINAL {'ADIABATIC' if self.opt.adiabatic else 'DIABATIC'} POPULATION: {final_obs['population']}")

        self._save(history, idx, nstate)

        return StateResult(
            population=final_obs["population"],
            kinetic_energy=final_obs["kinetic_energy"],
            potential_energy=final_obs["potential_energy"],
            total_energy=final_obs["total_energy"],
            position=final_obs["position"],
            momentum=final_obs["momentum"],
            norm=final_obs["norm"]
        )

    def _obs(self, wfn: Wavefunction, log: bool, pop_decay: np.ndarray, time: float = 0, wfn_0: Optional[Wavefunction] = None) -> dict:
        results = {}

        wfn_obs = wfn.to_adiabatic(self.grid, self.ham, time) if self.opt.adiabatic else wfn

        def should(f): return log or getattr(self.opt.write, f)

        if (should("autocorrelation") or self.opt.write.spectrum) and wfn_0:
            results["autocorrelation"] = wfn_0.overlap(wfn)
        if should("kinetic_energy") or should("total_energy"):
            results["kinetic_energy"] = wfn.ke(self.grid, self.ham)
        if should("momentum"):
            results["momentum"] = wfn.momentum(self.grid)
        if should("norm") or self.opt.absorbing_potential:
            results["norm"] = wfn.norm()
        if should("population"):
            results["population"] = wfn_obs.population(pop_decay)
        if should("position"):
            results["position"] = wfn.position(self.grid)
        if should("potential_energy") or should("total_energy"):
            results["potential_energy"] = wfn.pe(self.grid, self.ham, time)
        if should("total_energy"):
            results["total_energy"] = results["kinetic_energy"] + results["potential_energy"]
        if should("wavefunction"):
            results["wavefunction"] = wfn_obs.data.copy()

        return results

    def _pack_wfns(self, wfns: np.ndarray) -> np.ndarray:
        length, *_, nstate = wfns.shape

        grid = [r.ravel() for r in self.grid.position]
        wfns = wfns.reshape(length, -1, nstate).swapaxes(0, 1)

        return np.column_stack((*grid, wfns.view(float).reshape(len(wfns), -1)))

    def _get_spectrum(self, acf: np.ndarray) -> np.ndarray:
        times = np.arange(-len(acf) + 1, len(acf)) * self.opt.time_step
        tau = -np.log(1e-4) / (len(acf) * self.opt.time_step)**2

        acf = np.concatenate((acf[1:][::-1], np.conj(acf))) * np.exp(-tau * times**2)
        acf = np.pad(acf, 2 * [10 * len(acf)], mode="constant")

        omega = 2 * np.pi * np.fft.fftfreq(len(acf), self.opt.time_step)
        spectrum = np.fft.fftshift(np.fft.fft(np.fft.ifftshift(acf)) * self.opt.time_step)

        return np.column_stack((np.fft.fftshift(omega), np.real(spectrum)))

    def _head(self, idx: int, wfn: Wavefunction):
        ic, mode = self.opt.initial_conditions, "IMAGINARY" if self.opt.imaginary else "REAL"

        with np.printoptions(formatter={"float": "{:10.4f}".format}, suppress=True):
            print(f"\nSTATE {idx} INITIAL GAMMA: {np.array(ic.gamma, float)}\n")

        print(f"STATE {idx} {mode} TIME PROPAGATION")

        p_w, s_w = 11 * wfn.ndim + 1, 11 * wfn.nstate + 1
        pop_label = "ADIABATIC POP." if self.opt.adiabatic else "POPULATION"

        print(
            f"{'ITER':>7} {'KIN (Eh)':>12} {'POT (Eh)':>12} {'TOT (Eh)':>12} "
            f"{'POS (a0)':>{p_w}} {'MOM (hb/a0)':>{p_w}} "
            f"{pop_label:>{s_w}} {'NORM':>9} TIME"
        )

    def _log_step(self, i: int, obs: dict, duration: datetime.timedelta):
        with np.printoptions(formatter={"float": "{:10.4f}".format}, suppress=True):
            print(
                f"{i:7d} {obs['kinetic_energy']:12.6f} "
                f"{obs['potential_energy']:12.6f} {obs['total_energy']:12.6f} "
                f"{obs['position']} {obs['momentum']} {obs['population']} "
                f"{obs['norm']:1.3e} {duration}"
            )

    def _save(self, history: dict, idx: int, nstate: int):
        if not history: return
        
        times = np.arange(len(next(iter(history.values())))) * self.opt.time_step

        for field, path in self.opt.write:

            if not path or field not in history or not history[field]: continue

            field_arr, base, ext = np.array(history[field]), *path.rsplit(".", 1)

            if "autocorrelation" in field:
                data = np.column_stack((times, np.real(field_arr), np.imag(field_arr)))
            if "spectrum" in field:
                data = field_arr[0]
            elif "wavefunction" in field:
                data = self._pack_wfns(field_arr)
            else:
                data = np.column_stack((times, field_arr))
                
            fname = path if nstate == 1 else f"{base}_STATE-{idx:02}.{ext}"
            np.savetxt(fname, data, header=f"{data.shape[0]} {data.shape[1]}", comments="", fmt="%20.14f")
