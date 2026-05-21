from dataclasses import dataclass, field
from types import SimpleNamespace
from typing import Any

import numpy as np

from .grid import Grid


@dataclass(kw_only=True)
class History:
    autocorrelation: list[complex | None] = field(default_factory=list)
    kinetic_energy: list[float | None] = field(default_factory=list)
    momentum: list[np.ndarray | None] = field(default_factory=list)
    norm: list[float | None] = field(default_factory=list)
    population: list[np.ndarray | None] = field(default_factory=list)
    position: list[np.ndarray | None] = field(default_factory=list)
    potential_energy: list[float | None] = field(default_factory=list)
    total_energy: list[float | None] = field(default_factory=list)
    wavefunction: list[np.ndarray | None] = field(default_factory=list)

    @property
    def latest(self) -> SimpleNamespace:
        return SimpleNamespace(**{k: v[-1] if v else None for k, v in self.__dict__.items() if isinstance(v, list)})

    def get(self, name: str, grid: Grid, dt: float) -> np.ndarray:
        times = lambda: [np.arange(len(getattr(self, name))) * dt]

        grid_map = {
            "wavefunction": lambda: [r.ravel() for r in grid.pos]
        }

        data = np.array(getattr(self, name))

        if name == "wavefunction":
            data = np.moveaxis(data, 0, -2).reshape(-1, data.shape[0] * data.shape[-1])

        data = data.reshape(data.shape[0], -1)

        if np.iscomplexobj(data):
            data = np.ascontiguousarray(data).view(np.real(data).dtype)

        return np.column_stack((*grid_map.get(name, times)(), data))

    def record(self, **kwargs: Any) -> None:
        for k in kwargs.keys() & vars(self).keys(): getattr(self, k).append(kwargs[k])
