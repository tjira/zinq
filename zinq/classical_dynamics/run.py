import datetime
import time

from ..backend import np
from ..potential import Potential
from .options import Options
from .results import RunResult, TrajectoryResult
from .velocity_verlet import VelocityVerlet


def run(options_dict: dict):
    opt = Options(**options_dict)
    return Runner(opt).run()


class Runner:
    opt: Options
    pot: Potential

    def __init__(self, opt: Options):
        self.opt = opt
        self.potential = opt.potential.create()

    def run(self):
        pass
