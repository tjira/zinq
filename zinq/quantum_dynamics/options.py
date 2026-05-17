from typing import Optional, Union

from pydantic import BaseModel, Field

from ..potential import HarmonicOptions, TimeLinearOptions, TullyFirstOptions


class AbsorbingPotentialOptions(BaseModel):
    limits: list[list[float]]
    exponent: float = 0.001
    stop_norm: float = 1e-12
    track_population: bool = True


class GridOptions(BaseModel):
    limits: list[list[float]]
    npoint: int


class ImaginaryOptions(BaseModel):
    nstate: int = 1


class InitialConditionsOptions(BaseModel):
    adiabatic: bool = False
    gamma: list[float]
    momentum: list[float]
    position: list[float]
    state: int = 0


class WriteOptions(BaseModel):
    autocorrelation: Optional[str] = None
    final_wavefunction: Optional[str] = None
    kinetic_energy: Optional[str] = None
    momentum: Optional[str] = None
    norm: Optional[str] = None
    population: Optional[str] = None
    position: Optional[str] = None
    potential_energy: Optional[str] = None
    spectrum: Optional[str] = None
    total_energy: Optional[str] = None
    wavefunction: Optional[str] = None


class Options(BaseModel):
    absorbing_potential: Optional[AbsorbingPotentialOptions] = None
    adiabatic: bool = False
    grid: GridOptions
    imaginary: Optional[ImaginaryOptions] = None
    initial_conditions: InitialConditionsOptions
    iterations: Optional[int] = None
    log_interval: int = 1
    mass: float = 1
    time_step: float
    write: WriteOptions = WriteOptions()
    potential: Union[
        HarmonicOptions,
        TimeLinearOptions,
        TullyFirstOptions,
    ] = Field(discriminator="type")
