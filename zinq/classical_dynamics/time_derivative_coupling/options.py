from typing import Literal

from pydantic import BaseModel, ConfigDict

from .hammes_schiffer_tully import HammesSchifferTully


class HammesSchifferTullyOptions(BaseModel):
    model_config = ConfigDict(extra="forbid")

    type: Literal["hammes_schiffer_tully"]

    def create(self):
        return HammesSchifferTully()
