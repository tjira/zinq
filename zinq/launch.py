import argparse, json, os, pathlib, platform, subprocess, sys, tempfile

def executeInput(inp):
    with tempfile.NamedTemporaryFile(mode="w") as tmp:
        json.dump(inp, tmp); tmp.flush(); sys.exit(subprocess.call([getBinaryPath(), tmp.name]))

def getBinaryPath():
    ARCH, OS = platform.uname().machine.lower(), platform.uname().system.lower()

    if OS   == "darwin": OS   =   "macos"
    if ARCH ==  "arm64": ARCH = "aarch64"
    if ARCH ==  "amd64": ARCH =  "x86_64"

    binary_path = pathlib.Path(__file__).with_name("bin") / ("zinq-" + ARCH + "-" + OS + (".exe" if OS == "windows" else "")) 

    if OS != "windows": os.chmod(binary_path, 0o755)

    return binary_path

def getInputTemplate(name):
    return {
        "zinq" : [
            {
                "name": name,
                "options" : {}
            }
        ]
    }

def main():
    sys.exit(subprocess.call([getBinaryPath(), *sys.argv[1:]]))

def hf():
    parser = argparse.ArgumentParser(
        prog="Zinq Hartree-Fock Module", description="Wrapper for the Hartree-Fock method using the Zinq package.",
        formatter_class=lambda prog: argparse.HelpFormatter(prog, max_help_position=128),
        add_help=False, allow_abbrev=False
    )

    parser.add_argument("-h", "--help", action="help", default=argparse.SUPPRESS, help="This help message.")

    parser.add_argument("-b", "--basis", type=str, help="Basis set to use.", default="sto-3g")
    parser.add_argument("-m", "--molecule", type=str, help="Molecule specification in XYZ format.", default="molecule.xyz")

    args = parser.parse_args()

    if not os.path.isfile(args.molecule): raise Exception(f"ERROR: MOLECULE FILE '{args.molecule}' DOES NOT EXIST")

    inp = getInputTemplate("hartree_fock")

    inp["zinq"][0]["options"] = {
        "system" : args.molecule,
        "basis" : args.basis
    }

    executeInput(inp)

def prime():
    parser = argparse.ArgumentParser(
        prog="Zinq Prime Generation Module", description="Wrapper for the prime number generation using the Zinq package.",
        formatter_class=lambda prog: argparse.HelpFormatter(prog, max_help_position=128),
        add_help=False, allow_abbrev=False
    )

    parser.add_argument("-h", "--help", action="help", default=argparse.SUPPRESS, help="This help message.")

    parser.add_argument("-c", "--count", type=str, help="Number of primes to generate.", default=10)
    parser.add_argument("-l", "--log", type=str, help="Logging interval.", default=1)
    parser.add_argument("-o", "--output", type=str, help="Output file. Splitting interval can be put after a semicolon.")
    parser.add_argument("-s", "--start", type=int, help="Starting number.", default=2)

    parser.add_argument('--mersenne', action=argparse.BooleanOptionalAction)

    args = parser.parse_args()

    inp = getInputTemplate("prime_numbers")

    inp["zinq"][0]["options"] = {
        "count" : args.count,
        "filter" : "mersenne" if args.mersenne else "all",
        "log_interval" : args.log,
        "start" : args.start,
        "output" : {
            "interval" : int(args.output.split(":")[1]) if ":" in args.output else None,
            "path" : args.output.split(":")[0]
        } if args.output else None
    }

    executeInput(inp)
