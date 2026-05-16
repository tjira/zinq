from typing import Optional

from pydantic import BaseModel

from .landau_zener import LandauZenerOptions
from .surface_hopping import SurfaceHopping


class SurfaceHoppingOptions(BaseModel):
    landau_zener: Optional[LandauZenerOptions] = None

    def create(self) -> Optional[SurfaceHopping]:
        try: return next(opt.create() for _, opt in self if opt is not None)
        except StopIteration: return None
