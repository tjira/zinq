from typing import Optional

from ..backend import np
from ..potential import Potential
from .options import AbsorbingPotentialOptions


class Hamiltonian:
    pot: Potential
    mass: float
    cap: Optional[AbsorbingPotentialOptions]

    def __init__(self, pot: Potential, mass: float, cap: Optional[AbsorbingPotentialOptions] = None):
        self.pot, self.mass, self.cap = pot, mass, cap

        assert self.mass > 0, f"MASS MUST BE POSITIVE, GOT {self.mass}"

    def eval_cap(self, position: list[np.ndarray]) -> np.ndarray:
        if self.cap is None: return np.zeros_like(position[0])

        stiffness = self.cap.exponent
        total_penalty = np.zeros_like(position[0])

        for x, limits in zip(position, self.cap.limits):
            lower_penetration = np.maximum(0, limits[0] - x)
            upper_penetration = np.maximum(0, x - limits[1])
            
            step_penalty = np.exp(stiffness * (lower_penetration + upper_penetration)) - 1
            
            total_penalty += step_penalty

        return total_penalty
