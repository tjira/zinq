import argparse
import datetime
import os
import time

from zinq import __version__
from zinq.dispatcher import process_file


def main():
    parser = argparse.ArgumentParser(
        add_help=False,
        allow_abbrev=False,
        description="Lightweight Python engine for electronic structure and quantum dynamics.",
        formatter_class=lambda prog: argparse.HelpFormatter(prog, max_help_position=128),
        prog="zinq",
    )

    parser.add_argument("-h", "--help", action="help", help="This help message.")

    parser.add_argument(
        "inputs",
        nargs="?",
        default="input.json",
        help="Input files to process."
    )

    args = parser.parse_args()

    print(
        f"PYTHON: {os.sys.version.split()[0]}, ZINQ: {__version__}, "
        f"TIMESTAMP: {datetime.datetime.now().isoformat()}"
    )

    start_time = time.time()

    for input_file in args.inputs.split():
        process_file(input_file)

    duration = datetime.timedelta(seconds=time.time() - start_time)

    print(f"\nTOTAL EXECUTION TIME: {duration}")


if __name__ == "__main__":
    main()
