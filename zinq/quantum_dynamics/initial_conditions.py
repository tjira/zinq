from dataclasses import dataclass

import numpy as np

from .options import InitialConditionsConfig


@dataclass(frozen=True, kw_only=True)
class InitialConditions:
    adia: bool = False
    gamma: np.ndarray
    mom: np.ndarray
    pos: np.ndarray
    state: int

    @classmethod
    def from_options(cls, opt: InitialConditionsConfig):
        return cls(
            pos=np.array(opt.position),
            mom=np.array(opt.momentum),
            gamma=np.array(opt.gamma),
            state=opt.state,
            adia=opt.adiabatic,
        )
