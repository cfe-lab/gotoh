#!/usr/bin/env python3
""" The FastAPI web server
"""

import os
# import typing
# import pydantic.error_wrappers
from fastapi import FastAPI
# BackgroundTasks, Depends, Request, Form, WebSocket, HTTPException

from alignment.gotoh.about import __title__, __version__
from alignment.webconfig import WebServerConfig
import alignment.webtypes as webtypes
from alignment.align_it import AlignIt

API_VERSION = 1.0
API_DATE = "2022-02-02"


VersionResponse = webtypes.VersionResponse
AlignItInput = webtypes.AlignItInput
AlignItResult = webtypes.AlignItResult
AlignItResponse = webtypes.AlignItResponse
                    
app = FastAPI()

WS_CONFIG = WebServerConfig()


@app.get("/version",
         tags=["version"],
         response_model=VersionResponse)
def get_version() -> VersionResponse:
    """Return some version information about the fastapi server."""
    build_date = os.getenv("BUILD_DATE","Not found")
    return VersionResponse(api_version=API_VERSION,
                           api_date=API_DATE,
                           module_name = __title__,
                           module_version = __version__,
                           container_build_date = build_date)


@app.put("/align-it",
         tags=["alignment"],
         response_model=AlignItResponse)
def align_it(alinput: AlignItInput) -> AlignItResponse:
    """
    Do a sequence alignment and return the result.
    NOTE: Defaults to HyPhy pair scores - mismatch=5 mismatch=-4
          To use the old (Ruby) scores (1, -1), set: use_rb_match_scores=True
    """
    standard = alinput.seqa
    seq = alinput.seqb
    gap_ini = alinput.gap_init
    gap_ext = alinput.gap_extend
    use_terminal_gap_penalty = alinput.use_terminal_gap_penalty
    use_rb_match_scores = alinput.use_rb_match_scores 

    aligner = AlignIt()
    
    [sa, sb, score, sts] = aligner.align_it(standard, seq, 
                                            gap_ini, gap_ext, 
                                            use_terminal_gap_penalty,
                                            use_rb_match_scores)

    return AlignItResponse(seqa=sa,
                           seqb=sb,
                           align_score=score,
                           align_status=sts)


@app.put("/align-it-aa",
         tags=["alignment-aa"],
         response_model=AlignItResponse)
def align_it_aa(alinput: AlignItInput) -> AlignItResponse:
    """
    Do a sequence alignment and return the result.
    NOTE: Defaults to empirical score matrix based on 25% divergent HIV sequences
          See (Nickle et al., 2007)
          To use the old (Ruby) scores (4, 2), set: use_rb_match_scores=True
    """
    standard = alinput.seqa
    seq = alinput.seqb
    gap_ini = alinput.gap_init
    gap_ext = alinput.gap_extend
    use_terminal_gap_penalty = alinput.use_terminal_gap_penalty
    use_rb_match_scores = alinput.use_rb_match_scores

    aligner = AlignIt()

    [sa, sb, score, sts] = aligner.align_it_aa(standard, seq, 
                                               gap_ini, gap_ext, 
                                               use_terminal_gap_penalty,
                                               use_rb_match_scores)

    return AlignItResponse(seqa=sa,
                           seqb=sb,
                           align_score=score,
                           align_status=sts)
