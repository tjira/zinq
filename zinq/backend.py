import numpy as np
import scipy.linalg as linalg


def enable_cupy():
    import cupy
    import cupyx.scipy.linalg

    globals()["np"] = cupy
    globals()["linalg"] = cupyx.scipy.linalg
