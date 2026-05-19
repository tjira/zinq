from dataclasses import dataclass
from typing import Optional

from ..potential import Potential


@dataclass
class Hamiltonian:
    pot: Potential
    cap: Optional[Potential]
    mass: float
