from ..potential import Potential


class Hamiltonian:
    pot: Potential
    mass: float

    def __init__(self, pot: Potential, mass: float):
        self.pot, self.mass = pot, mass

        assert self.mass > 0, f"MASS MUST BE POSITIVE, GOT {self.mass}"
