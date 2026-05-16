import datetime
import time

from ..backend import np
from ..potential import Potential
from .options import Options
from .ensemble import Ensemble
from .results import RunResult, TrajectoryResult
from .velocity_verlet import VelocityVerlet


def run(options_dict: dict):
    opt = Options(**options_dict)
    return Runner(opt).run()


class Runner:
    opt: Options
    potential: Potential

    def __init__(self, opt: Options):
        self.opt, self.potential = opt, opt.potential.create()

    def run(self):
        ensemble = Ensemble(
            np.array(self.opt.initial_conditions.position),
            np.array(self.opt.initial_conditions.momentum),
            np.array(self.opt.initial_conditions.gamma),
            self.opt.trajectories,
            self.opt.initial_conditions.state
        )

        verlet = VelocityVerlet(self.potential, self.opt.mass, self.opt.time_step, 1e-8)

        for i in range(self.opt.iterations + 1):
            time_current = i * self.opt.time_step

            if i > 0: verlet.step(ensemble, time_current - self.opt.time_step)

            log = i == 0 or i == self.opt.iterations or (i % self.opt.log_interval == 0)

            kinetic_energy = ensemble.ke()
            potential_energy = ensemble.pe(self.potential)

            if log:
                print(kinetic_energy, potential_energy, kinetic_energy + potential_energy)
