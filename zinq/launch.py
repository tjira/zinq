import os, pathlib, platform, subprocess, sys

def main():
    ARCH, OS = platform.uname().machine.lower(), platform.uname().system.lower()

    if OS   == "darwin": OS   =   "macos"
    if ARCH ==  "arm64": ARCH = "aarch64"
    if ARCH ==  "amd64": ARCH =  "x86_64"

    binary_path = pathlib.Path(__file__).with_name("bin") / ("zinq-" + ARCH + "-" + OS + (".exe" if OS == "windows" else "")) 

    if OS != "windows": os.chmod(binary_path, 0o755)

    exit_code = subprocess.call([binary_path, *sys.argv[1:]])

    sys.exit(exit_code)
