from dataclasses import dataclass
import numpy as np

@dataclass(frozen=True, kw_only=True)
class TullyFirst:
    A: float = 0.01
    B: float = 1.6
    C: float = 0.005
    D: float = 1

    def evaluateDiabatic(self, r: np.ndarray):
        V00 = np.sign(r) * self.A * (1 - np.exp(-self.B * np.abs(r)))
        V01 = self.C * np.exp(-self.D * r**2)
        V22 = np.sign(r) * self.A * (np.exp(-self.B * np.abs(r)) - 1)

        return np.array([[V00, V01], [V01, V22]])
