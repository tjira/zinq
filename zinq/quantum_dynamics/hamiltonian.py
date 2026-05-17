from typing import Optional

from ..backend import np
from ..potential import Potential
from .options import ComplexAbsorbingPotentialOptions


class Hamiltonian:
    pot: Potential
    mass: float
    cap: Optional[ComplexAbsorbingPotentialOptions]

    def __init__(self, pot: Potential, mass: float, cap: Optional[ComplexAbsorbingPotentialOptions] = None):
        self.pot, self.mass, self.cap = pot, mass, cap

        assert self.mass > 0, f"MASS MUST BE POSITIVE, GOT {self.mass}"

    def eval_cap(self, position: list[np.ndarray]) -> np.ndarray:
        if self.cap is None: return np.zeros_like(position[0])

        exp, zipped = self.cap.exponent, zip(position, self.cap.limits)

        return sum((np.exp(exp * (np.maximum(0, l[0] - x) + np.maximum(0, x - l[1]))) - 1 for x, l in zipped), start=np.zeros_like(position[0]))
