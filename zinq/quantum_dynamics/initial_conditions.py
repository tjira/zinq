from dataclasses import dataclass
import numpy as np

from typing_extensions import TypedDict


@dataclass(frozen=True, kw_only=True)
class InitialConditions:
    Config = TypedDict("Config", {
        "position": list[float],
        "momentum": list[float],
        "gamma": list[float],
        "state": int,
        "adiabatic": bool,
    })

    position: list[float]
    momentum: list[float]
    gamma: list[float]
    state: int = 0
    adiabatic: bool = False

    def __post_init__(self):
        object.__setattr__(self, "pos", np.array(self.position, dtype=float))
        object.__setattr__(self, "mom", np.array(self.momentum, dtype=float))
        object.__setattr__(self, "gamma", np.array(self.gamma, dtype=float))
