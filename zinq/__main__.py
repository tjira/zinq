import argparse
import datetime
import os
import time

from zinq import __version__


def main():
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

    parser.add_argument("inputs", nargs="?", default="input.json", help="Input files.")

    args = parser.parse_args()

    if args.cupy:
        from . import backend
        backend.enable_cupy()

    from .dispatcher import process_file

    print(
        f"PYTHON: {os.sys.version.split()[0]}, ZINQ: {__version__}, "
        f"TIMESTAMP: {datetime.datetime.now().isoformat()}, "
        f"BACKEND: {'CUPY' if args.cupy else 'NUMPY'}"
    )

    if args.profile:
        import cProfile, pstats, subprocess
        profiler = cProfile.Profile()
        profiler.enable()

    start_time = time.time()

    for input_file in args.inputs.split():
        process_file(input_file)

    duration = datetime.timedelta(seconds=time.time() - start_time)

    print(f"\nTOTAL EXECUTION TIME: {duration}")

    if args.profile:
        profiler.disable()

        stats_file, dot_file, pdf_file = "profile.stats", "profile.dot", "profile.pdf"

        gprof2dot_cmd = [
            os.sys.executable, "-m", "gprof2dot",
            "-f", "pstats",
            "-n", "1.0",
            "-e", "0.5",
            stats_file,
            "-o", dot_file
        ]

        dot_cmd = [
            "dot",
            "-Tpdf",
            dot_file, 
            "-o", pdf_file
        ]

        profiler.dump_stats(stats_file)
        
        subprocess.run(gprof2dot_cmd, check=True)
        subprocess.run(dot_cmd, check=True)


if __name__ == "__main__":
    main()
