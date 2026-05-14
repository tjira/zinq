import numpy as np
from .split_operator import SplitOperator
from .wavefunction import Wavefunction
from ..potential import TullyFirst

def run(options: dict):
    pot = TullyFirst()
    print(pot)
