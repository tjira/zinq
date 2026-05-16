from typing import Optional

from pydantic import BaseModel

from ..potential import PotentialOptions


class InitialConditionsOptions(BaseModel):
    position: list[float]
    momentum: list[float]
    gamma: list[float]
    state: int = 0


class LogIntervalOptions(BaseModel):
    iteration: int = 1
    trajectory: int = 1


class WriteOptions(BaseModel):
    trajectory: Optional[str] = None
    energy: Optional[str] = None


class Options(BaseModel):
    initial_conditions: InitialConditionsOptions
    potential: PotentialOptions
    iterations: Optional[int] = None
    trajectories: Optional[int] = None
    time_step: float
    mass: float = 1
    log_interval: LogIntervalOptions = LogIntervalOptions()
    write: WriteOptions = WriteOptions()
