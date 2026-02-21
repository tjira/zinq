import os, platform, setuptools, setuptools.command.build_py, shutil, subprocess

ARCH, OS = platform.uname().machine.lower().replace("arm64", "aarch64").replace("amd64", "x86_64") , platform.uname().system.lower().replace("darwin", "macos")

def version(fallback="0.0.0"):

    if not shutil.which("git"): return fallback

    try:
        result = subprocess.run(["git", "describe", "--tags"], capture_output=True, text=True)

        if result.returncode != 0: return fallback

        return result.stdout.strip().lstrip("v").replace("-", "+", 1).replace("-", ".")

    except Exception: return fallback

class Build(setuptools.command.build_py.build_py):
    def run(self):

        environment = {**os.environ, **({"OS" : "Windows_NT"} if OS == "windows" else {})}

        subprocess.run(["make", "CROSS=1"], check=True, env=environment)

        super().run()

binaries = [f"zig-out/{ARCH}-{OS}/" + binary + (".exe" if OS == "windows" else "") for binary in ["zinq"]]

setuptools.setup(
    version = version(),
    packages = setuptools.find_packages(),
    has_ext_modules=lambda: True,
    data_files = [("Scripts" if OS == "windows" else "bin", binaries)],
    cmdclass = {"build_py": Build}
)
