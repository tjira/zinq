from typing import Optional

from pydantic import BaseModel

from ..potential import PotentialOptions
from .grid import GridOptions
from .hamiltonian import ComplexAbsorbingPotentialOptions


class ImaginaryOptions(BaseModel):
    nstate: int = 1


class InitialConditionsOptions(BaseModel):
    position: list[float]
    momentum: list[float]
    gamma: list[float]
    state: int = 0
    adiabatic: bool = False


class WriteOptions(BaseModel):
    autocorrelation: Optional[str] = None
    kinetic_energy: Optional[str] = None
    momentum: Optional[str] = None
    norm: Optional[str] = None
    population: Optional[str] = None
    position: Optional[str] = None
    potential_energy: Optional[str] = None
    total_energy: Optional[str] = None
    wavefunction: Optional[str] = None
    spectrum: Optional[str] = None
    final_wavefunction: Optional[str] = None


class Options(BaseModel):
    grid: GridOptions
    initial_conditions: InitialConditionsOptions
    potential: PotentialOptions
    iterations: Optional[int] = None
    time_step: float
    mass: float = 1
    log_interval: int = 1
    adiabatic: bool = False
    imaginary: Optional[ImaginaryOptions] = None
    absorbing_potential: Optional[ComplexAbsorbingPotentialOptions] = None
    write: WriteOptions = WriteOptions()
