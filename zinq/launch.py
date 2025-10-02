import os, pathlib, platform, sys

def main():
    ARCH, OS = platform.uname().machine.lower(), platform.uname().system.lower()

    if OS   == "darwin": OS   =   "macos"
    if ARCH ==  "arm64": ARCH = "aarch64"

    binary_path = pathlib.Path(__file__).with_name("bin") / ("zinq-" + ARCH + "-" + OS + (".exe" if OS == "windows" else "")) 

    if OS != "windows": os.chmod(binary_path, 0o755)

    os.execv(binary_path, [binary_path, *sys.argv[1:]])
