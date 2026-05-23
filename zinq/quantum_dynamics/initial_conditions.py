"""Module defining the initial conditions for quantum dynamics simulations."""

from dataclasses import dataclass

import numpy as np

from .options import InitialConditionsConfig


@dataclass(frozen=True, kw_only=True)
class InitialConditions:
    """
    Initial conditions for a quantum dynamics simulation wavepacket.

    Attributes
    ----------
    adia : bool
        Whether the initial state is represented in the adiabatic basis.
    gamma : np.ndarray
        Width parameters of the initial Gaussian wavepacket.
    mom : np.ndarray
        Initial momentum components of the wavepacket.
    pos : np.ndarray
        Initial center position of the wavepacket.
    state : int
        Initial electronic state index.

    """

    adia: bool = False
    gamma: np.ndarray
    mom: np.ndarray
    pos: np.ndarray
    state: int

    @classmethod
    def from_options(cls, opt: InitialConditionsConfig) -> "InitialConditions":
        """
        Create an InitialConditions instance from a configuration object.

        Parameters
        ----------
        opt : InitialConditionsConfig
            Configuration options for initial conditions.

        Returns
        -------
        InitialConditions
            The initialized initial conditions object.

        """
        return cls(
            pos=np.array(opt.position),
            mom=np.array(opt.momentum),
            gamma=np.array(opt.gamma),
            state=opt.state,
            adia=opt.adiabatic,
        )
