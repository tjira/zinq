import os
import platform
import subprocess
import tarfile
import urllib.request
import zipfile
from pathlib import Path

from setuptools import Extension, setup
from setuptools.command.build_ext import build_ext

ZIG_VERSION = "0.16.0"


class ZigExtension(Extension):
    def __init__(self, name: str, source: str):
        super().__init__(name, sources=[source])


class ZigBuild(build_ext):
    def build_extension(self, ext: Extension) -> None:
        if not isinstance(ext, ZigExtension):
            return super().build_extension(ext)

        dest_path = os.path.abspath(self.get_ext_fullpath(ext.name))
        flags = ["-dynamic", f"-femit-bin={dest_path}", "--name", ext.name.split(".")[-1], "-O", "ReleaseFast"]

        os.makedirs(os.path.dirname(dest_path), exist_ok=True)
        subprocess.check_call([self._ensure_zig(), "build-lib", *flags, ext.sources[0]])

    def _download_and_extract_zig(self, local_dir: Path, system: str, machine: str, ext: str) -> None:
        url = f"https://ziglang.org/download/{ZIG_VERSION}/zig-{machine}-{system}-{ZIG_VERSION}.{ext}"
        archive_path = local_dir / f"zig.{ext}"

        local_dir.mkdir(parents=True, exist_ok=True)
        urllib.request.urlretrieve(url, archive_path)

        if ext == "tar.xz":
            with tarfile.open(archive_path, "r:xz") as tar:
                tar.extractall(local_dir)

        if ext == "zip":
            with zipfile.ZipFile(archive_path, "r") as zip_ref:
                zip_ref.extractall(local_dir)

        if archive_path.exists():
            archive_path.unlink()

    def _get_zig_platform_details(self) -> tuple[str, str, str]:
        system_name = platform.system().lower()
        machine_name = platform.machine().lower()

        system = {"darwin": "macos"}.get(system_name, system_name)
        machine = {"amd64": "x86_64", "arm64": "aarch64"}.get(machine_name, machine_name)

        return system, machine, "zip" if system == "windows" else "tar.xz"

    def _ensure_zig(self) -> str:
        local_dir = Path("build/compiler")

        if path := self._find_zig_exe(local_dir):
            return path

        system, machine, ext = self._get_zig_platform_details()

        self._download_and_extract_zig(local_dir, system, machine, ext)

        if path := self._find_zig_exe(local_dir):
            return path

        raise FileNotFoundError("ZIG EXECUTABLE COULD NOT BE FOUND AFTER DOWNLOAD AND EXTRACTION")

    def _find_zig_exe(self, local_dir: Path) -> str | None:
        if not local_dir.exists(): return None
    
        exe = "zig.exe" if platform.system() == "Windows" else "zig"

        return next((str(p) for p in local_dir.rglob(exe) if os.access(p, os.X_OK)), None)


if __name__ == "__main__":
    os.makedirs("zinq/lib", exist_ok=True)

    modules = [
        ZigExtension("zinq.lib.native", "zinq/zig/native.zig"),
    ]

    setup(ext_modules=modules, cmdclass={"build_ext": ZigBuild})
