import json
import os

from . import classical_dynamics, quantum_dynamics


def process_file(file_path):
    if not os.path.exists(file_path):
        raise FileNotFoundError(f"INPUT FILE '{file_path}' NOT FOUND")

    with open(file_path) as f:
        input_json = json.load(f)

    print(f"\nFILE TO PROCESS: {file_path}\n\n{json.dumps(input_json, indent=2)}")

    for item in input_json["zinq"]:
        match item["name"]:
            case "quantum_dynamics":
                result = quantum_dynamics.run(item["options"])
            case "classical_dynamics":
                result = classical_dynamics.run(item["options"])
