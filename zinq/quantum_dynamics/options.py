from typing import Literal

from pydantic import BaseModel, ConfigDict

from ..potential.potential import AnyPotential


class AbsorberConfig(BaseModel):
    model_config = ConfigDict(extra="forbid")

    limits: list[tuple[float, float]]
    exponent: float = 0.001


class GridConfig(BaseModel):
    model_config = ConfigDict(extra="forbid")

    limits: list[tuple[float, float]]
    npoint: int


class HamiltonianConfig(BaseModel):
    model_config = ConfigDict(extra="forbid")

    potential: AnyPotential
    mass: float
    absorber: AbsorberConfig | None = None


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
    imaginary: ImaginaryConfig | None = None
    initial_conditions: InitialConditionsConfig
    iterations: int | None = None
    log_interval: int = 1
    stop_norm: float = 1e-8
    time_step: float
    write: WriteConfig | None = None
