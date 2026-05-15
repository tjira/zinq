from ..potential import Potential


class Hamiltonian:
    potential: Potential
    mass: float

    def __init__(self, potential: Potential, mass: float):
        self.potential, self.mass = potential, mass

        assert self.mass > 0, f"MASS MUST BE POSITIVE, GOT {self.mass}"
