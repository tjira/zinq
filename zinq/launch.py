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

def eigh():
    parser = argparse.ArgumentParser(
        prog="Zinq Eigenproblem Solver Module", description="Wrapper for the eigenproblem solver using the Zinq package.",
        formatter_class=lambda prog: argparse.HelpFormatter(prog, max_help_position=128),
        add_help=False, allow_abbrev=False
    )

    parser.add_argument("-h", "--help", action="help", default=argparse.SUPPRESS, help="This help message.")

    parser.add_argument("-c", "--eigenvectors", type=str, help="Path to the file where to write the eigenvectors.")
    parser.add_argument("-i", "--input", type=str, help="Path to the input matrix.", default="A.mat")
    parser.add_argument("-j", "--eigenvalues", type=str, help="Path to the file where to write the eigenvalues.")

    parser.add_argument("--print", action=argparse.BooleanOptionalAction, default=False, help="Print the input matrix, eigenvalues, and eigenvectors to the console.")

    args = parser.parse_args()

    inp = getInputTemplate("eigenproblem_solver")

    inp["zinq"][0]["options"] = {
        "matrix_file" : args.input,
        "print" : {
            "input_matrix" : args.print,
            "eigenvalues" : args.print,
            "eigenvectors" : args.print
        },
        "write" : {
            "eigenvalues" : args.eigenvalues if args.eigenvalues else None,
            "eigenvectors" : args.eigenvectors if args.eigenvectors else None
        },
        "hermitian" : True,
        "real" : True
    }

    executeInput(inp)

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

    inp = getInputTemplate("hartree_fock")

    inp["zinq"][0]["options"] = {
        "system" : args.molecule,
        "basis" : args.basis
    }

    executeInput(inp)

def molint():
    parser = argparse.ArgumentParser(
        prog="Zinq Molecular Integrals Module", description="Wrapper for the calculation of molecular integrals using the Zinq package.",
        formatter_class=lambda prog: argparse.HelpFormatter(prog, max_help_position=128),
        add_help=False, allow_abbrev=False
    )

    parser.add_argument("-h", "--help", action="help", default=argparse.SUPPRESS, help="This help message.")

    parser.add_argument("-b", "--basis", type=str, help="Basis set to use.", default="sto-3g")
    parser.add_argument("-j", "--coulomb", type=str, help="Molecule specification in XYZ format.", default="")
    parser.add_argument("-m", "--molecule", type=str, help="Molecule specification in XYZ format.", default="molecule.xyz")
    parser.add_argument("-s", "--overlap", type=str, help="Calculation of overlap integrals.", default="")
    parser.add_argument("-t", "--kinetic", type=str, help="Calculation of kinetic energy integrals.", default="")
    parser.add_argument("-v", "--nuclear", type=str, help="Calculation of nuclear attraction integrals.", default="")

    args = parser.parse_args()

    inp = getInputTemplate("molecular_integrals")

    inp["zinq"][0]["options"] = {
        "system" : args.molecule,
        "basis" : args.basis,
        "write" : {
            "coulomb" : args.coulomb if args.coulomb else None,
            "overlap" : args.overlap if args.overlap else None,
            "kinetic" : args.kinetic if args.kinetic else None,
            "nuclear" : args.nuclear if args.nuclear else None
        }
    }

    executeInput(inp)

def mp2():
    parser = argparse.ArgumentParser(
        prog="Zinq MP2 Module", description="Wrapper for the MP2 method using the Zinq package.",
        formatter_class=lambda prog: argparse.HelpFormatter(prog, max_help_position=128),
        add_help=False, allow_abbrev=False
    )

    parser.add_argument("-h", "--help", action="help", default=argparse.SUPPRESS, help="This help message.")

    parser.add_argument("-b", "--basis", type=str, help="Basis set to use.", default="sto-3g")
    parser.add_argument("-m", "--molecule", type=str, help="Molecule specification in XYZ format.", default="molecule.xyz")

    args = parser.parse_args()

    inp = getInputTemplate("moller_plesset")

    inp["zinq"][0]["options"] = {
        "hartree_fock" : {
            "system" : args.molecule,
            "basis" : args.basis
        },
        "order" : 2
    }

    executeInput(inp)

def primecheck():
    parser = argparse.ArgumentParser(
        prog="Zinq Prime Checking Module", description="Wrapper for the prime number check using the Zinq package.",
        formatter_class=lambda prog: argparse.HelpFormatter(prog, max_help_position=128),
        add_help=False, allow_abbrev=False
    )

    parser.add_argument("-h", "--help", action="help", default=argparse.SUPPRESS, help="This help message.")

    parser.add_argument("-b", "--bits", type=int, help="Number of bits to use.", default=64)
    parser.add_argument("-n", "--number", type=str, help="Number to check.", required=True)

    parser.add_argument("--mersenne", action=argparse.BooleanOptionalAction, default=False, help="Check only Mersenne primes.")

    args = parser.parse_args()

    inp = getInputTemplate("prime_numbers")

    inp["zinq"][0]["options"] = {
        "mode" : {
            "check" : {
                "filter" : "mersenne" if args.mersenne else "all",
                "number" : args.number
            }
        },
        "bits" : args.bits if args.bits else None
    }

    executeInput(inp)

def primefact():
    parser = argparse.ArgumentParser(
        prog="Zinq Prime Factorization Module", description="Wrapper for the prime number factorization using the Zinq package.",
        formatter_class=lambda prog: argparse.HelpFormatter(prog, max_help_position=128),
        add_help=False, allow_abbrev=False
    )

    parser.add_argument("-h", "--help", action="help", default=argparse.SUPPRESS, help="This help message.")

    parser.add_argument("-b", "--bits", type=int, help="Number of bits to use.", default=64)
    parser.add_argument("-n", "--number", type=str, help="Number to factorize.", required=True)

    args = parser.parse_args()

    inp = getInputTemplate("prime_numbers")

    inp["zinq"][0]["options"] = {
        "mode" : {
            "factorize" : {
                "number" : args.number
            }
        },
        "bits" : args.bits if args.bits else None
    }

    executeInput(inp)

def primegen():
    parser = argparse.ArgumentParser(
        prog="Zinq Prime Generation Module", description="Wrapper for the prime number generation using the Zinq package.",
        formatter_class=lambda prog: argparse.HelpFormatter(prog, max_help_position=128),
        add_help=False, allow_abbrev=False
    )

    parser.add_argument("-h", "--help", action="help", default=argparse.SUPPRESS, help="This help message.")

    parser.add_argument("-b", "--bits", type=int, help="Number of bits to use.", default=64)
    parser.add_argument("-c", "--count", type=int, help="Number of primes to generate.", default=10)
    parser.add_argument("-l", "--log", type=str, help="Logging interval.", default=1)
    parser.add_argument("-o", "--output", type=str, help="Output file. Splitting interval can be put after a semicolon.")
    parser.add_argument("-s", "--start", type=str, help="Starting number.", default="2")

    parser.add_argument('--mersenne', action=argparse.BooleanOptionalAction)

    args = parser.parse_args()

    inp = getInputTemplate("prime_numbers")

    inp["zinq"][0]["options"] = {
        "mode" : {
            "generate" : {
                "count" : args.count,
                "filter" : "mersenne" if args.mersenne else "all",
                "log_interval" : args.log,
                "start" : args.start,
                "output" : {
                    "interval" : int(args.output.split(":")[1]) if ":" in args.output else None,
                    "path" : args.output.split(":")[0]
                } if args.output else None
            }
        },
        "bits" : args.bits if args.bits else None
    }

    executeInput(inp)
