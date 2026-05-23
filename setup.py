"""Setup script for compiling Zig extensions and building the zinq package."""

import subprocess
import sys
from pathlib import Path

from setuptools import Extension, setup
from setuptools.command.build_ext import build_ext


class ZigExtension(Extension):
    """
    Representation of a Zig C-compatible extension module.

    Parameters
    ----------
    name : str
        The full module name.
    source : str
        The path to the Zig source file.

    """

    def __init__(self, name: str, source: str) -> None:
        """
        Initialize the ZigExtension.

        Parameters
        ----------
        name : str
            Module name.
        source : str
            Zig source file path.

        """
        super().__init__(name, sources=[source])


class ZigBuild(build_ext):
    """Custom build extension command to compile Zig source code into shared libraries."""

    def build_extension(self, ext: Extension) -> None:
        """
        Build a ZigExtension compile target.

        Parameters
        ----------
        ext : Extension
            The extension target to build.

        """
        if not isinstance(ext, ZigExtension):
            super().build_extension(ext)
            return

        dest_path = Path(self.get_ext_fullpath(ext.name)).resolve()
        flags = [
            "-dynamic",
            f"-femit-bin={dest_path}",
            "--name",
            ext.name.split(".")[-1],
            "-O",
            "ReleaseFast",
        ]

        dest_path.parent.mkdir(parents=True, exist_ok=True)
        subprocess.check_call(  # noqa: S603
            [sys.executable, "-m", "ziglang", "build-lib", *flags, ext.sources[0]]
        )


if __name__ == "__main__":
    Path("zinq/lib").mkdir(parents=True, exist_ok=True)

    modules = [
        ZigExtension("zinq.lib.native", "zinq/zig/native.zig"),
    ]

    setup(ext_modules=modules, cmdclass={"build_ext": ZigBuild})
