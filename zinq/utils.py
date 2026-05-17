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

    try:
        import cupy
        package_versions["CUPY"] = cupy.__version__
    except ImportError:
        pass

    return package_versions


def print_startup_header(backend_name=None):
    package_versions = get_versions()
    
    header = f"PYTHON: {sys.version.split()[0]}, ZINQ: {__version__}"
    header += f", TIMESTAMP: {datetime.datetime.now().isoformat()}"

    if backend_name: header += f", BACKEND: {backend_name}\n"
    
    print(header)

    for pkg, ver in package_versions.items():
        print(f"{pkg}: {ver}")
