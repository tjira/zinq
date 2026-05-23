"""Module for loading Zig dynamic libraries."""

import ctypes
from pathlib import Path


def load_library(lib_name: str) -> ctypes.CDLL:
    """
    Load a compiled Zig shared library by name.

    Parameters
    ----------
    lib_name : str
        The name of the Zig library to load.

    Returns
    -------
    ctypes.CDLL
        The loaded dynamic library object.

    """
    lib_dir = Path(__file__).parent / "../lib"
    lib_path = next(iter(lib_dir.glob(f"{lib_name}.*")))
    return ctypes.CDLL(str(lib_path))
