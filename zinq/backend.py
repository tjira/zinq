import numpy as np
import scipy.linalg as linalg


def enable_cupy():
    global np, linalg

    import cupy
    import cupyx.scipy.linalg

    np, linalg = cupy, cupyx.scipy.linalg
