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
    seqa: str = ""
    seqb: str = ""
    gap_init: int = 6
    gap_extend: int = 1
    use_terminal_gap_penalty: bool = False
    use_rb_match_scores: bool = False


class AlignItResult(str, Enum):
    """Return status of the align_it function"""
    ok = "ok"
    illegal_char = "illegal_char"
    internal_error = "internal_error"


class AlignItResponse(BaseModel):
    """The response of the align_it function"""

    class Config:
        use_enum_values = True

    seqa: str = ""
    seqb: str = ""
    align_score: int = 0
    align_status: AlignItResult
