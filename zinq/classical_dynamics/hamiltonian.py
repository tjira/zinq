from ..potential import Potential


class Hamiltonian:
    pot: Potential
    mass: float

    def __init__(self, pot: Potential, mass: float):
        assert mass > 0, "MASS MUST BE POSITIVE"

        self.pot, self.mass = pot, mass
