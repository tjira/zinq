from typing import Optional

from pydantic import BaseModel

from ..backend import np
from ..potential import Potential


class ComplexAbsorbingPotentialOptions(BaseModel):
    limits: list[list[float]]
    exponent: float = 0.001
    stop_norm: float = 1e-12
    track_population: bool = True


class Hamiltonian:
    potential: Potential
    mass: float
    cap: Optional[ComplexAbsorbingPotentialOptions]

    def __init__(self, potential: Potential, mass: float, cap: Optional[ComplexAbsorbingPotentialOptions] = None):
        self.potential, self.mass, self.cap = potential, mass, cap

        assert self.mass > 0, f"MASS MUST BE POSITIVE, GOT {self.mass}"

    def eval_cap(self, position: list[np.ndarray]) -> np.ndarray:
        if self.cap is None: return np.zeros_like(position[0])

        exp, zipped = self.cap.exponent, zip(position, self.cap.limits)

        return sum((np.exp(exp * (np.maximum(0, l[0] - x) + np.maximum(0, x - l[1]))) - 1 for x, l in zipped), start=np.zeros_like(position[0]))
