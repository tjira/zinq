import argparse
import cProfile
import datetime
import shutil
import subprocess
import sys
import time

from zinq import __version__

from . import backend
from .utils import print_startup_header


def main():
    profiler = None
    parser = argparse.ArgumentParser(
        add_help=False,
        allow_abbrev=False,
        description="Lightweight Python engine for electronic structure and quantum dynamics.",
        formatter_class=lambda prog: argparse.HelpFormatter(prog, max_help_position=128),
        prog="zinq",
    )

    parser.add_argument("-h", "--help", action="help", help="This help message.")
    parser.add_argument("--cupy", action="store_true", help="Use Cupy instead of Numpy.")
    parser.add_argument("--profile", action="store_true", help="Profile the execution and generate a PDF graph.")

    parser.add_argument("inputs", nargs="*", default=["input.json"], help="Input files.")

    args = parser.parse_args()

    if args.cupy: backend.enable_cupy()

    from .dispatcher import process_file

    print_startup_header(backend_name="CUPY" if args.cupy else "NUMPY")

    if args.profile:
        profiler = cProfile.Profile()
        profiler.enable()

    start_time = time.time()

    print(f"\nINPUT FILES TO PROCESS: {args.inputs}")

    for input_file in args.inputs: process_file(input_file)

    duration = datetime.timedelta(seconds=time.time() - start_time)

    print(f"\nTOTAL EXECUTION TIME: {duration}")

    if profiler:
        profiler.disable()

        stats_file, dot_file, pdf_file = "profile.stats", "profile.dot", "profile.svg"

        gprof2dot_cmd = [
            sys.executable, "-m", "gprof2dot",
            "-f", "pstats",
            "-n", "1.0",
            "-e", "0.5",
            stats_file,
            "-o", dot_file
        ]

        dot_cmd = [
            "dot",
            "-Tsvg",
            dot_file, 
            "-o", pdf_file
        ]

        profiler.dump_stats(stats_file)

        if not shutil.which("gprof2dot"): raise RuntimeError("EXECUTABLE 'gprof2dot' NOT FOUND")
        
        subprocess.run(gprof2dot_cmd, check=True)

        if not shutil.which("dot"): raise RuntimeError("EXECUTABLE 'dot' NOT FOUND")

        subprocess.run(dot_cmd, check=True)


if __name__ == "__main__":
    main()
