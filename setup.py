import os, platform, setuptools, setuptools.command.build_py, shutil, subprocess, setuptools.command.bdist_wheel

def archos():

    replacements = {
        "arm64": "aarch64",
        "amd64": "x86_64",
        "darwin": "macos"
    }

    ARCH, OS = platform.uname().machine.lower(), platform.uname().system.lower()

    PLATFORM = os.environ.get("PLATFORM", f"{replacements.get(OS, OS)}_{replacements.get(ARCH, ARCH)}")

    return PLATFORM

def version(fallback="0.0.0"):

    if not shutil.which("git"): return fallback

    try:
        result = subprocess.run(["git", "describe", "--tags"], capture_output=True, text=True)

        if result.returncode != 0: return fallback

        return result.stdout.strip().lstrip("v").replace("-", "+", 1).replace("-", ".")

    except Exception: return fallback

PLATFORM = archos(); OS, ARCH = PLATFORM.split("_", 1)

class Bdist(setuptools.command.bdist_wheel.bdist_wheel):
    def get_tag(self):

        replacements = {
            "linux" : "manylinux_2_5",
            "macos" : "macosx_10_0",
            "windows" : "win"
        }

        plat = f"{replacements.get(OS, OS)}_{replacements.get(ARCH, ARCH)}"

        return "py3", "none", plat.replace("win_aarch64", "win_arm64").replace("win_x86_64", "win_amd64")

class Build(setuptools.command.build_py.build_py):
    def run(self):

        environment = {**os.environ, **({"OS" : "Windows_NT"} if os.name == "nt" else {})}

        subprocess.run(["make", "CROSS=1"], check=True, env=environment)

        super().run()

binaries = [f"zig-out/{ARCH}-{OS}/" + binary + (".exe" if OS == "windows" else "") for binary in ["zinq"]]

setuptools.setup(
    version = version(),
    packages = setuptools.find_packages(),
    has_ext_modules=lambda: True,
    data_files = [("Scripts" if OS == "windows" else "bin", binaries)],
    cmdclass = {"bdist_wheel": Bdist, "build_py": Build}
)
