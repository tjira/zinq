from typing import Optional

from pydantic import BaseModel

from ..backend import np
from ..potential import Potential
from .options import ComplexAbsorbingPotentialOptions


class Hamiltonian:
    potential: Potential
    mass: float
    cap: Optional[ComplexAbsorbingPotentialOptions]

    def __init__(self, potential: Potential, mass: float, cap: Optional[ComplexAbsorbingPotentialOptions] = None):
        self.potential, self.mass, self.cap = potential, mass, cap

        assert self.mass > 0, f"MASS MUST BE POSITIVE, GOT {self.mass}"

    def eval_cap(self, position: list[np.ndarray]) -> np.ndarray:
        if self.cap is None: return 0

        return sum(np.exp(self.cap.exponent * (np.maximum(0, l[0] - x) + np.maximum(0, x - l[1]))) - 1 for x, l in zip(position, self.cap.limits))
