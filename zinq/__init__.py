import os, platform, sys


def _get_binary_path():
    suffix = ".exe" if platform.system() == "Windows" else ""

    return os.path.join(os.path.dirname(__file__), "bin", "zinq" + suffix)


def main():
    os.execvp(_get_binary_path(), [_get_binary_path()] + sys.argv[1:])
