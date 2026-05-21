from typing import Literal

from pydantic import BaseModel, ConfigDict

from ..potential.potential import AnyPotential


class GridConfig(BaseModel):
    model_config = ConfigDict(extra="forbid")

    limits: list[tuple[float, float]]
    npoint: int


class HamiltonianConfig(BaseModel):
    model_config = ConfigDict(extra="forbid")

    potential: AnyPotential
    mass: float


class ImaginaryConfig(BaseModel):
    model_config = ConfigDict(extra="forbid")

    nstate: int


class InitialConditionsConfig(BaseModel):
    model_config = ConfigDict(extra="forbid")

    position: list[float]
    momentum: list[float]
    gamma: list[float]
    state: int = 0
    adiabatic: bool = False


class WriteConfig(BaseModel):
    model_config = ConfigDict(extra="forbid")

    autocorrelation: str | None = None
    kinetic_energy: str | None = None
    momentum: str | None = None
    norm: str | None = None
    population: str | None = None
    position: str | None = None
    potential_energy: str | None = None
    spectrum: str | None = None
    total_energy: str | None = None
    wavefunction: str | None = None


class Options(BaseModel):
    model_config = ConfigDict(extra="forbid")
    type: Literal["quantum_dynamics"]

    adiabatic: bool = False
    grid: GridConfig
    hamiltonian: HamiltonianConfig
    initial_conditions: InitialConditionsConfig
    imaginary: ImaginaryConfig | None = None
    iterations: int | None = None
    log_interval: int = 1
    time_step: float
    write: WriteConfig | None = None
