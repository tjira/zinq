import os, setuptools, setuptools.command.install, shutil, subprocess, sysconfig

from setuptools.command.build_py import build_py as _build_py

class Build(_build_py):
    def run(self):
        super().run()

        subprocess.run(["make", "CROSS=1"], check=True)

        os.makedirs(os.path.join("zinq", "bin"), exist_ok=True)

        for directory in os.listdir("zig-out"):

            source = os.path.join("zig-out", directory, "zinq.exe" if "win" in directory else "zinq")
            dest = os.path.join("zinq", "bin", "zinq-" + directory + (".exe" if "win" in directory else ""))

            shutil.copyfile(source, dest)

setuptools.setup(
    name = "zinq",
    version = open("VERSION").read(),
    packages = setuptools.find_packages(),
    long_description = open("README.md").read(),
    long_description_content_type = "text/markdown",
    author = "tjira",
    url = "https://github.com/tjira/zinq",
    cmdclass = {
        "build_py": Build
    },
    entry_points = {
        "console_scripts": ["zinq = zinq.launch:main"]
    },
    package_data = {
        "zinq": ["bin/*"]
    }
)
