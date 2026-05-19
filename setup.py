import os
import subprocess

from setuptools import Extension, setup
from setuptools.command.build_ext import build_ext


class ZigExtension(Extension):
    def __init__(self, name, source):
        super().__init__(name, sources=[source])


class ZigBuild(build_ext):
    def build_extension(self, ext):
        if not isinstance(ext, ZigExtension):
            return super().build_extension(ext)

        dest_path = os.path.abspath(self.get_ext_fullpath(ext.name))

        flags = [
            "-dynamic", f"-femit-bin={dest_path}", "--name", ext.name.split(".")[-1], "-O", "ReleaseFast"
        ]

        os.makedirs(os.path.dirname(dest_path), exist_ok=True)

        subprocess.check_call(["zig", "build-lib", *flags, ext.sources[0]])


if __name__ == "__main__":
    os.makedirs("zinq/lib", exist_ok=True)
        
    modules = [
        ZigExtension("zinq.lib.native", "zinq/zig/native.zig"),
    ]

    setup(ext_modules=modules, cmdclass={"build_ext": ZigBuild})
