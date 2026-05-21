from dataclasses import dataclass, field

import numpy as np

from .grid import Grid
from .observables import Observables


@dataclass(kw_only=True)
class History:
    obs: list[Observables] = field(default_factory=lambda: [])

    def append(self, obs: Observables) -> None:
        self.obs.append(obs)

    def get(self, name: str, grid: Grid, dt: float) -> np.ndarray:
        times = lambda: [np.arange(len(self.obs)) * dt]

        grid_map = {
            "wavefunction": lambda: [r.ravel() for r in grid.pos]
        }

        data = np.array([getattr(obs, name) for obs in self.obs])

        if name == "wavefunction":
            data = np.moveaxis(data, 0, -2).reshape(-1, data.shape[0] * data.shape[-1])
            
        data = data.reshape(data.shape[0], -1)

        if np.iscomplexobj(data):
            data = np.ascontiguousarray(data).view(np.real(data).dtype)

        return np.column_stack((*grid_map.get(name, times)(), data))
