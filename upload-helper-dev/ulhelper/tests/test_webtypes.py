"""Test the webtypes module"""


import ulhelper.webtypes as webtypes
import pytest
from pydantic.error_wrappers import ValidationError


class TestWebTypes:
    """Test the webtypes"""

    def test_versionresponse01(self) -> None:
        """VersionResponse() without arguments
        raises a ValidationError"""
        with pytest.raises(ValidationError):
            webtypes.VersionResponse()

    def test_versionresponse02(self) -> None:
        """We can instantiate a VersionResponse
        instance with correct arguments."""
        vr = webtypes.VersionResponse(api_version='1',
                                      api_date='1900-08-04')
        assert vr is not None, "vr is None"
