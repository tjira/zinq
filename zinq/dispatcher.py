import json
import os
from typing import Annotated, List, Union

from pydantic import BaseModel, ConfigDict, Field

from . import quantum_dynamics


class Options(BaseModel):
    model_config = ConfigDict(extra="forbid")

    zinq: List[Annotated[Union[
        quantum_dynamics.Options
    ], Field(discriminator="type")]]


def process_file(file_path):
    if not os.path.exists(file_path):
        raise FileNotFoundError(f"INPUT FILE '{file_path}' NOT FOUND")

    with open(file_path) as f:
        input_json = json.load(f)

    print(f"\nFILE TO PROCESS: {file_path}\n\n{json.dumps(input_json, indent=2)}")

    for item in Options.model_validate(input_json).zinq:

        match item.type:
            case "quantum_dynamics":
                quantum_dynamics.run(item) 
