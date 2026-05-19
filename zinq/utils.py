import datetime
import sys

from zinq import __version__


def get_versions():
    import matplotlib
    import numpy
    import pydantic
    import scipy
    import sympy

    package_versions = {
        "NUMPY": numpy.__version__,
        "SCIPY": scipy.__version__,
        "PYDANTIC": pydantic.__version__,
        "MATPLOTLIB": matplotlib.__version__,
        "SYMPY": sympy.__version__,
    }

    return package_versions


def print_startup_header():
    package_versions = get_versions()
    
    header = f"PYTHON: {sys.version.split()[0]}, ZINQ: {__version__}"
    header += f", TIMESTAMP: {datetime.datetime.now().isoformat()}\n"
    
    print(header)

    for pkg, ver in package_versions.items():
        print(f"{pkg}: {ver}")
