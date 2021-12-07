#!/usr/bin/env python3
""" The FastAPI web server
"""

# import typing
# import pydantic.error_wrappers
from fastapi import FastAPI
# BackgroundTasks, Depends, Request, Form, WebSocket, HTTPException

from ulhelper.webconfig import WebServerConfig
import ulhelper.webtypes as webtypes


API_VERSION = 1.0
API_DATE = "2021-12-06"


VersionResponse = webtypes.VersionResponse
AlignItInput = webtypes.AlignItInput
AlignItResult = webtypes.AlignItResult
AlignItResponse = webtypes.AlignItResponse

app = FastAPI()

WS_CONFIG = WebServerConfig()


@app.get("/hello", tags=["test"])
def hullo():
    """A simple hello world test. Does not require DB access."""
    return {"Hello": "World"}


@app.get("/version",
         tags=["version"],
         response_model=VersionResponse)
def get_version() -> VersionResponse:
    """Return some version information about the fastapi server.
    Does not require DB access."""
    return VersionResponse(api_version=API_VERSION,
                           api_date=API_DATE)


@app.put("/align-it",
         tags=["alignment"],
         response_model=AlignItResponse)
def align_it(ainput: AlignItInput) -> AlignItResponse:
    """Do a sequence alignment and return the result.
    NOTE: if this takes too long, we will have to implement
    running this in a background job.
    """
    return AlignItResponse(seqb='bananas',
                           align_status=AlignItResult.internal_error)
