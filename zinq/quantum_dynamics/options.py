from typing import Optional

from pydantic import BaseModel

from ..potential import HarmonicOptions, Potential, TimeLinearOptions, TullyFirstOptions


class ComplexAbsorbingPotentialOptions(BaseModel):
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
    position: list[float]
    momentum: list[float]
    gamma: list[float]
    state: int = 0
    adiabatic: bool = False


class PotentialOptions(BaseModel):
    harmonic: Optional[HarmonicOptions] = None
    time_linear: Optional[TimeLinearOptions] = None
    tully_1: Optional[TullyFirstOptions] = None

    def create(self) -> Potential:
        try: return next(pot.create() for _, pot in self if pot is not None)
        except StopIteration: raise ValueError("NO POTENTIAL SPECIFIED")


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
