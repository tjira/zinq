import os, platform, sys, subprocess


def _get_binary_path():
    suffix = ".exe" if platform.system() == "Windows" else ""

    return os.path.join(os.path.dirname(__file__), "bin", "zinq" + suffix)


def main():
    sys.exit(subprocess.call([_get_binary_path()] + sys.argv[1:]))
