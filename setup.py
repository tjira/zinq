import os, setuptools, setuptools.command.build_py, shutil, subprocess

class Build(setuptools.command.build_py.build_py):
    def run(self):
        subprocess.run(["make", "CROSS=1"], check=True, env={**os.environ, **({"OS" : "Windows_NT"} if os.name == "nt" else {})})

        os.makedirs(os.path.join("zinq", "bin"), exist_ok=True)

        for directory in os.listdir("zig-out"):

            suffix = ".exe" if "windows" in directory else ""

            source, dest = os.path.join("zig-out", directory, "zinq" + suffix), os.path.join("zinq", "bin", "zinq-" + directory + suffix)

            shutil.copyfile(source, dest)

        super().run()

setuptools.setup(
    version = open("VERSION").read().strip(),
    packages = setuptools.find_packages(),
    cmdclass = {
        "build_py": Build
    }
)
