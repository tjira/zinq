from ...backend import np
from .time_derivative_coupling import TimeDerivativeCoupling


class HammesSchifferTully(TimeDerivativeCoupling):
    def evaluate(self, U_old: np.ndarray, U_new: np.ndarray, dt: float) -> np.ndarray:
        overlap = np.einsum("nki,nkj->nij", U_old, U_new)
        return (overlap - np.moveaxis(overlap, [1, 2], [2, 1])) / (2 * dt)
