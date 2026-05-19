import ctypes
from pathlib import Path


def load_library(lib_name: str):
    return ctypes.CDLL(str(list((Path(__file__).parent / "../lib").glob(f"{lib_name}.*"))[0]))
