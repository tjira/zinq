import argparse
import contextlib
import cProfile
import datetime
import shutil
import subprocess
import sys
import time

from zinq import __version__

from .utils import load_library, print_startup_header


def main():
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

    from .dispatcher import process_file

    print_startup_header()

    with cProfile.Profile() if args.profile else contextlib.nullcontext() as profiler:
        start_time = time.time()

        print(f"\nINPUT FILES TO PROCESS: {args.inputs}")

        for input_file in args.inputs: process_file(input_file)

        duration = datetime.timedelta(seconds=time.time() - start_time)

        print(f"\nTOTAL EXECUTION TIME: {duration}")

    if args.profile:
        stats_file, dot_file, svg_file = "profile.stats", "profile.dot", "profile.svg"

        profiler.dump_stats(stats_file) # type: ignore

        if not shutil.which("gprof2dot"): raise RuntimeError("EXECUTABLE 'gprof2dot' NOT FOUND")
        if not shutil.which("dot"): raise RuntimeError("EXECUTABLE 'dot' NOT FOUND")

        subprocess.run([
            sys.executable,
            "-m", "gprof2dot",
            "-f", "pstats",
            "-n", "1.0",
            "-e", "0.5",
            stats_file,
            "-o", dot_file
        ], check=True)

        subprocess.run(["dot", "-Tsvg", dot_file, "-o", svg_file], check=True)

        print(f"ZIG ADD(10, 32) = {load_library('native').add(10, 32)}")


if __name__ == "__main__":
    main()
