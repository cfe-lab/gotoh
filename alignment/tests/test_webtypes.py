"""Test the webtypes module"""


import alignment.webtypes as webtypes
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

    def test_alignitinput01(self) -> None:
        """We can instantiate a AlignItInput
        instance with correct arguments."""
        aip = webtypes.AlignItInput(seqa="gattaca",
                                    seqb="gattaa",
                                    gap_init=6,
                                    gap_extend=1,
                                    use_terminal_gap_penalty=False,
                                    use_rb_match_scores=False)
        assert aip is not None, "aip is None"

    def test_alignitresult01(self) -> None:
        """AlignItResult() without specifying value
        raises a TypeError"""
        with pytest.raises(TypeError):
            webtypes.AlignItResult()

    def test_alignitresponse01(self) -> None:
        """We can instantiate a AlignItResponse
        instance with correct arguments."""
        airp = webtypes.AlignItResponse(seqa="gattaca",
                                        seqb="gatta-a",
                                        align_score=6,
                                        align_status=webtypes.AlignItResult.ok)
        assert airp is not None, "airp is None"