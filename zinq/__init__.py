import os, platform, sys


def _get_binary_path():
    return os.path.join(os.path.dirname(__file__), "bin", "zinq")


def main():
    os.execvp(_get_binary_path(), [_get_binary_path()] + sys.argv[1:])
