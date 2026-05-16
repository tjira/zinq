import datetime
import time

from ..backend import np
from ..potential import Potential
from .options import Options
from .ensemble import Ensemble
from .results import RunResult
from .velocity_verlet import VelocityVerlet


def run(options_dict: dict):
    opt = Options(**options_dict)
    return Runner(opt).run()


class Runner:
    opt: Options
    potential: Potential

    def __init__(self, opt: Options):
        self.opt, self.potential = opt, opt.potential.create()

    def run(self) -> RunResult:
        ensemble = Ensemble(
            np.array(self.opt.initial_conditions.position),
            np.array(self.opt.initial_conditions.momentum),
            np.array(self.opt.initial_conditions.gamma),
            self.opt.trajectories,
            self.opt.initial_conditions.state,
            self.opt.seed
        )

        verlet = VelocityVerlet(self.potential, self.opt.mass, self.opt.time_step, 1e-8)

        history = {f: [] for f, p in self.opt.write if p is not None}
        
        self._head(ensemble)

        for i in range(self.opt.iterations + 1):
            start, current_time = time.time(), i * self.opt.time_step

            if i > 0: verlet.step(ensemble, current_time - self.opt.time_step)

            log = i == 0 or i == self.opt.iterations or (i % self.opt.log_interval == 0)

            obs = self._obs(ensemble, log, current_time)

            for f, v in obs.items():
                if f in history: history[f].append(v)

            if log: self._log_step(i, obs, datetime.timedelta(seconds=time.time() - start))

        self._save(history)

        final = self._obs(ensemble, True, self.opt.iterations * self.opt.time_step)

        print(f"\nFINAL ADIABATIC POPULATION: {final['population']}")

        return RunResult(**{k: final[k] for k in RunResult.__annotations__})

    def _head(self, ensemble: Ensemble):
        ic = self.opt.initial_conditions
        
        with np.printoptions(formatter={"float": "{:10.4f}".format}, suppress=True):
            print(f"\nINITIAL GAMMA: {np.array(ic.gamma, float)}\n")
            
        print(f"CLASSICAL ENSEMBLE TIME PROPAGATION")

        p_w = 11 * ensemble.ndim + 1
        s_w = 11 * self.potential.nstate + 1

        print(
            f"{'ITER':>7} {'KIN (Eh)':>12} {'POT (Eh)':>12} {'TOT (Eh)':>12} "
            f"{'POS (a0)':>{p_w}} {'MOM (hb/a0)':>{p_w}} "
            f"{'POPULATION':>{s_w}} TIME"
        )

    def _log_step(self, i: int, obs: dict, duration: datetime.timedelta):
        with np.printoptions(formatter={"float": "{:10.4f}".format}, suppress=True):
            print(
                f"{i:7d} {obs['kinetic_energy']:12.6f} {obs['potential_energy']:12.6f} {obs['total_energy']:12.6f} "
                f"{obs['position']} {obs['momentum']} {obs['population']} {duration}"
            )

    def _obs(self, ensemble: Ensemble, log: bool, time: float = 0) -> dict:
        results, write = {}, self.opt.write

        def has(f): return log or getattr(write, f)

        if has("population"): results["population"] = ensemble.population(self.potential.nstate)
        if has("position"): results["position"] = ensemble.position()
        if has("momentum"): results["momentum"] = ensemble.momentum()
        if has("kinetic_energy") or has("total_energy"): results["kinetic_energy"] = ensemble.ke(self.opt.mass)
        if has("potential_energy") or has("total_energy"): results["potential_energy"] = ensemble.pe(self.potential, time)
        if has("total_energy"): results["total_energy"] = results["kinetic_energy"] + results["potential_energy"]

        return results

    def _save(self, history: dict):
        times = np.arange(len(next(iter(history.values())))) * self.opt.time_step if history else []

        for field, path in self.opt.write:
            if not path or field not in history or not history[field]: continue

            data = np.column_stack((times, np.array(history[field])))
            np.savetxt(path, data, header=f"{data.shape[0]} {data.shape[1]}", comments="", fmt="%20.14f")
