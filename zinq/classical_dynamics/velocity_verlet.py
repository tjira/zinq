from ..backend import np
from ..potential import Potential


class VelocityVerlet:
    pot: Potential
    mass: float


    def __init__(self, potential: Potential, mass: float):
        self.potential = potential
        self.mass = mass

    def step(self):
        pass
