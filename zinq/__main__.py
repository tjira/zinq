"""Main entry point for the zinq engine CLI."""

import argparse
import contextlib
import cProfile
import datetime
import json
import shutil
import subprocess
import sys
import time
from pathlib import Path
from typing import Annotated

import matplotlib as mpl
import numpy as np
import pydantic
import scipy
from pydantic import BaseModel, ConfigDict, Field

from zinq import __version__

from . import quantum_dynamics


class Options(BaseModel):
    """
    Schema for parsing simulation input configurations.

    Attributes
    ----------
    zinq : list[Annotated[quantum_dynamics.Options, Field(discriminator="type")]]
        List of simulations to run.

    """

    model_config = ConfigDict(extra="forbid")

    zinq: list[Annotated[quantum_dynamics.Options, Field(discriminator="type")]]


def main() -> None:
    """Run the main command-line interface of the zinq engine."""
    parser = argparse.ArgumentParser(
        add_help=False,
        allow_abbrev=False,
        description="Lightweight Python engine for electronic structure and quantum dynamics.",
        formatter_class=lambda prog: argparse.HelpFormatter(prog, max_help_position=128),
        prog="zinq",
    )

    help_help = "Show this help message and exit."
    profile_help = "Profile the execution and generate a PDF graph."

    parser.add_argument("-h", "--help", action="help", help=help_help)
    parser.add_argument("--profile", action="store_true", help=profile_help)

    parser.add_argument("inputs", nargs="*", default=["input.json"], help="Input files.")

    args = parser.parse_args()

    _print_startup_header()

    with cProfile.Profile() if args.profile else contextlib.nullcontext() as profiler:
        start_time = time.time()

        print(f"\nINPUT FILES TO PROCESS: {args.inputs}")

        for input_file in args.inputs:
            _process_file(input_file)

        duration = datetime.timedelta(seconds=time.time() - start_time)

        print(f"\nTOTAL EXECUTION TIME: {duration}")

    if args.profile:
        stats_file, dot_file, svg_file = "profile.stats", "profile.dot", "profile.svg"

        profiler.dump_stats(stats_file) if profiler else None

        gprof2dot_path = shutil.which("gprof2dot")
        if not gprof2dot_path:
            msg = "EXECUTABLE 'gprof2dot' NOT FOUND"
            raise RuntimeError(msg)

        dot_path = shutil.which("dot")
        if not dot_path:
            msg = "EXECUTABLE 'dot' NOT FOUND"
            raise RuntimeError(msg)

        subprocess.run(
            [
                sys.executable,
                "-m",
                "gprof2dot",
                "-f",
                "pstats",
                "-n",
                "1.0",
                "-e",
                "0.5",
                stats_file,
                "-o",
                dot_file,
            ],
            check=True,
        )

        subprocess.run([dot_path, "-Tsvg", dot_file, "-o", svg_file], check=True)


def _get_versions() -> dict[str, str]:
    """Get versions of major scientific library dependencies."""
    return {
        "NUMPY": np.__version__,
        "SCIPY": scipy.__version__,
        "PYDANTIC": pydantic.__version__,
        "MATPLOTLIB": mpl.__version__,
    }


def _print_startup_header() -> None:
    """Print system configuration and library versions at startup."""
    package_versions = _get_versions()

    header = f"PYTHON: {sys.version.split()[0]}, ZINQ: {__version__}"
    header += f", TIMESTAMP: {datetime.datetime.now(datetime.timezone.utc).isoformat()}\n"

    print(header)

    for pkg, ver in package_versions.items():
        print(f"{pkg}: {ver}")


def _process_file(file_path: str) -> None:
    """Process a single JSON input file and run the configured simulations."""
    path = Path(file_path)
    if not path.exists():
        msg = f"INPUT FILE '{file_path}' NOT FOUND"
        raise FileNotFoundError(msg)

    with path.open() as f:
        input_json = json.load(f)

    print(f"\nFILE TO PROCESS: {file_path}\n\n{json.dumps(input_json, indent=2)}")

    for item in Options.model_validate(input_json).zinq:
        match item.type:
            case "quantum_dynamics":
                quantum_dynamics.run(item)


if __name__ == "__main__":
    main()
