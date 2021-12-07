"""Define classes and types for web services
"""

from pydantic import BaseModel
from enum import Enum


class VersionResponse(BaseModel):
    """The response to an API version request."""
    api_version: str
    api_date: str


class AlignItInput(BaseModel):
    """The input parameters to the align_it function.
    seqb is aligned to seqa."""
    seqa: str
    seqb: str
    gap_init: int
    gap_penalty: int


class AlignItResult(str, Enum):
    """Return status of the Align_it function"""
    ok = "ok"
    illegal_char = "illegal_char"
    internal_error = "internal_error"


class AlignItResponse(BaseModel):
    """The response of the align_it function"""

    class Config:
        use_enum_values = True

    seqb: str = ""
    align_status: AlignItResult
