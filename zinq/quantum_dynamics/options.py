"""Configuration options and Pydantic schemas for quantum dynamics simulations."""

from typing import Literal

from pydantic import BaseModel, ConfigDict

from zinq.potential.union import AnyPotential


class AbsorberConfig(BaseModel):
    """
    Configuration options for the boundary absorbing potential.

    Attributes
    ----------
    limits : list[tuple[float, float]]
        The boundary limits for the absorbing potential along each dimension.
    exponent : float
        The exponent determining the absorption strength.

    """

    model_config = ConfigDict(extra="forbid")

    limits: list[tuple[float, float]]
    exponent: float = 0.001


class GridConfig(BaseModel):
    """
    Configuration options for the spatial simulation grid.

    Attributes
    ----------
    limits : list[tuple[float, float]]
        The spatial limits of the grid along each dimension.
    npoint : int
        The number of grid points per dimension.

    """

    model_config = ConfigDict(extra="forbid")

    limits: list[tuple[float, float]]
    npoint: int


class HamiltonianConfig(BaseModel):
    """
    Configuration options for the Hamiltonian operator.

    Attributes
    ----------
    potential : AnyPotential
        The potential energy operator description.
    mass : float
        The particle mass in atomic units.
    absorber : AbsorberConfig | None
        Optional boundary absorber configuration.

    """

    model_config = ConfigDict(extra="forbid")

    potential: AnyPotential
    mass: float
    absorber: AbsorberConfig | None = None


class ImaginaryConfig(BaseModel):
    """
    Configuration options for imaginary time propagation.

    Attributes
    ----------
    nstate : int
        The number of electronic states to consider in imaginary time.

    """

    model_config = ConfigDict(extra="forbid")

    nstate: int


class InitialConditionsConfig(BaseModel):
    """
    Configuration options for initial wavepacket conditions.

    Attributes
    ----------
    position : list[float]
        The initial coordinate center of the wavepacket.
    momentum : list[float]
        The initial momentum components of the wavepacket.
    gamma : list[float]
        The initial wavepacket width parameter(s).
    state : int
        The initial electronic state index.
    adiabatic : bool
        Whether the state is initialized in the adiabatic representation.

    """

    model_config = ConfigDict(extra="forbid")

    position: list[float]
    momentum: list[float]
    gamma: list[float]
    state: int = 0
    adiabatic: bool = False


class WriteConfig(BaseModel):
    """
    Configuration options for output file writing targets.

    Attributes
    ----------
    autocorrelation : str | None
        Path to save autocorrelation function data.
    kinetic_energy : str | None
        Path to save kinetic energy data.
    momentum : str | None
        Path to save expectation value of momentum.
    norm : str | None
        Path to save wavepacket norm over time.
    population : str | None
        Path to save state populations.
    position : str | None
        Path to save coordinate expectation values.
    potential_energy : str | None
        Path to save potential energy over time.
    spectrum : str | None
        Path to save spectrum data.
    total_energy : str | None
        Path to save total energy.
    wavefunction : str | None
        Path to save wavepacket snapshots.

    """

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
    """
    Root configuration options for the entire quantum dynamics simulation.

    Attributes
    ----------
    type : Literal["quantum_dynamics"]
        Discriminator field identifying the type of simulation.
    adiabatic : bool
        Whether to perform the propagation in the adiabatic representation.
    grid : GridConfig
        The simulation grid configuration.
    hamiltonian : HamiltonianConfig
        The Hamiltonian operator configuration.
    imaginary : ImaginaryConfig | None
        Optional configuration for imaginary time propagation.
    initial_conditions : InitialConditionsConfig
        The initial wavepacket configuration.
    iterations : int | None
        The number of propagation steps.
    log_interval : int
        Interval at which to log status to console/file.
    stop_norm : float
        Norm threshold below which simulation stops.
    time_step : float
        Real or imaginary propagation time step.
    write : WriteConfig | None
        Configuration for files to output.

    """

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
