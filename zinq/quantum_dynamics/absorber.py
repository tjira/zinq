from dataclasses import dataclass

import numpy as np

from .options import AbsorberConfig


@dataclass(frozen=True, kw_only=True)
class Absorber:
    limits: np.ndarray
    exponent: float

    @classmethod
    def from_options(cls, opt: AbsorberConfig):
        return cls(
            limits=np.array(opt.limits),
            exponent=opt.exponent,
        )

    def eval(self, r: list[np.ndarray]) -> np.ndarray:
        total_penalty = np.zeros_like(r[0])

        for x, limits in zip(r, self.limits):
            lower_penetration = np.maximum(0, limits[0] - x)
            upper_penetration = np.maximum(0, x - limits[1])
            
            step_penalty = np.exp(self.exponent * (lower_penetration + upper_penetration)) - 1
            
            total_penalty += step_penalty

        return total_penalty
