"""Test the align_it module"""


import alignment.align_it as al
import pytest
from pydantic.error_wrappers import ValidationError


class TestAlignIt:
    """Test align_it"""

    def test_emptyseqs(self) -> None:
        """align_it() with empty seqs returns empty seqs
        raises a ValidationError"""
        aligner = al.AlignIt()
        [sa, sb, score, sts] = aligner.align_it(seqa="",
                                                seqb="", 
                                                gap_ini=1,
                                                gap_ext=1,
                                                use_terminal_gap_penalty=False,
                                                emulate_rb=False)
        assert sa == "", "empty seqa expected"
        assert sb == "", "empty seqb expected"
        assert score == 0, "expected score=0"
        assert sts == "ok", "alignment failed"

    def test_emptyseqs_aa(self) -> None:
        """align_it_aa() with empty seqs returns empty seqs
        raises a ValidationError"""
        aligner = al.AlignIt()
        [sa, sb, score, sts] = aligner.align_it_aa(seqa="",
                                                   seqb="", 
                                                   gap_ini=1,
                                                   gap_ext=1,
                                                   use_terminal_gap_penalty=False,
                                                   emulate_rb=False)
        assert sa == "", "empty seqa expected"
        assert sb == "", "empty seqb expected"
        assert score == 0, "expected score=0"
        assert sts == "ok", "alignment failed"

    def test_emptyseq_a(self) -> None:
        """align_it() with empty seq a returns gaps
        raises a ValidationError"""
        aligner = al.AlignIt()
        [sa, sb, score, sts] = aligner.align_it(seqa="",
                                                   seqb="gattaca", 
                                                   gap_ini=6,
                                                   gap_ext=1,
                                                   use_terminal_gap_penalty=False,
                                                   emulate_rb=False)
        assert sa == "-------", "seqa incorrect"
        assert sb == "gattaca", "seqb incorrect"
        assert score == 13, "expected score=13"
        assert sts == "ok", "alignment failed"

    def test_emptyseq_b(self) -> None:
        """align_it() with empty seq b returns gaps
        raises a ValidationError"""
        aligner = al.AlignIt()
        [sa, sb, score, sts] = aligner.align_it(seqa="gattaca",
                                                   seqb="", 
                                                   gap_ini=6,
                                                   gap_ext=1,
                                                   use_terminal_gap_penalty=False,
                                                   emulate_rb=False)
        assert sa == "gattaca", "seqa incorrect"
        assert sb == "-------", "seqb incorrect"
        assert score == 13, "expected score=13"
        assert sts == "ok", "alignment failed"

    def test_emptyseq_a_aa(self) -> None:
        """align_it_aa() with empty seq a returns gaps
        raises a ValidationError"""
        aligner = al.AlignIt()
        [sa, sb, score, sts] = aligner.align_it_aa(seqa="",
                                                   seqb="gattaca", 
                                                   gap_ini=6,
                                                   gap_ext=1,
                                                   use_terminal_gap_penalty=False,
                                                   emulate_rb=False)
        assert sa == "-------", "seqa incorrect"
        assert sb == "gattaca", "seqb incorrect"
        assert score == 13, "expected score=13"
        assert sts == "ok", "alignment failed"

    def test_emptyseq_b_aa(self) -> None:
        """align_it_aa() with empty seq b returns gaps
        raises a ValidationError"""
        aligner = al.AlignIt()
        [sa, sb, score, sts] = aligner.align_it_aa(seqa="gattaca",
                                                   seqb="", 
                                                   gap_ini=6,
                                                   gap_ext=1,
                                                   use_terminal_gap_penalty=False,
                                                   emulate_rb=False)
        assert sa == "gattaca", "seqa incorrect"
        assert sb == "-------", "seqb incorrect"
        assert score == 13, "expected score=13"
        assert sts == "ok", "alignment failed"

    def test_invalidchars(self) -> None:
        """align_it() with invalid chars returns an error
        raises a ValidationError"""
        aligner = al.AlignIt()
        [_, _, _, sts] = aligner.align_it(seqa="gattaza",
                                          seqb="gattaca", 
                                          gap_ini=6,
                                          gap_ext=1,
                                          use_terminal_gap_penalty=False,
                                          emulate_rb=False)
        assert sts == "illegal_char", "alignment did not fail with illegal_char"

    def test_invalidchars_aa(self) -> None:
        """align_it_aa() with invalid chars returns an error
        raises a ValidationError"""
        aligner = al.AlignIt()
        [_, _, _, sts] = aligner.align_it_aa(seqa="gatta$a",
                                             seqb="gattaca", 
                                             gap_ini=6,
                                             gap_ext=1,
                                             use_terminal_gap_penalty=False,
                                             emulate_rb=False)
        assert sts == "illegal_char", "alignment did not fail with illegal_char"

    def test_int_err(self) -> None:
        """align_it_aa() with invalid args returns an error
        raises a ValidationError"""
        aligner = al.AlignIt()
        [_, _, _, sts] = aligner.align_it(seqa="gattaca",
                                          seqb= 42.0, 
                                          gap_ini=1,
                                          gap_ext=1,
                                          use_terminal_gap_penalty=False,
                                          emulate_rb=False)
        assert sts == "internal_error", "alignment did not fail with internal_error"

    def test_int_err_aa(self) -> None:
        """align_it_aa() with invalid args returns an error
        raises a ValidationError"""
        aligner = al.AlignIt()
        [_, _, _, sts] = aligner.align_it_aa(seqa="gattaca",
                                             seqb= 42.0, 
                                             gap_ini=1,
                                             gap_ext=1,
                                             use_terminal_gap_penalty=False,
                                             emulate_rb=False)
        assert sts == "internal_error", "alignment did not fail with internal_error"

    def test_align_it(self) -> None:
        """align_it_aa() with test sequences
        raises a ValidationError"""
        aligner = al.AlignIt()
        [sa, sb, score, sts] = aligner.align_it(seqa="gattacagattaca",
                                                seqb="gattacttaca", 
                                                gap_ini=6,
                                                gap_ext=1,
                                                use_terminal_gap_penalty=False,
                                                emulate_rb=False)
        assert sa == "gattacagattaca", "seqa incorrect"
        assert sb == "gattac---ttaca", "seqb incorrect"
        assert score == 46, "expected score=76"
        assert sts == "ok", "alignment failed"

    def test_align_it_rb(self) -> None:
        """align_it_aa_rb() with test sequences, legacy (ruby) score matrix
        raises a ValidationError"""
        aligner = al.AlignIt()
        [sa, sb, score, sts] = aligner.align_it(seqa="gattacagattaca",
                                                seqb="gattacttaca", 
                                                gap_ini=6,
                                                gap_ext=1,
                                                use_terminal_gap_penalty=False,
                                                emulate_rb=True)
        assert sa == "gattacagattaca----", "seqa incorrect"
        assert sb == "-------gattacttaca", "seqb incorrect"
        assert score == 0, "expected score=0"
        assert sts == "ok", "alignment failed"

    def test_align_it_aa(self) -> None:
        """align_it_aa() with test sequences
        raises a ValidationError"""
        aligner = al.AlignIt()
        [sa, sb, score, sts] = aligner.align_it_aa(seqa="gattacagattaca",
                                                   seqb="gattacttaca", 
                                                   gap_ini=6,
                                                   gap_ext=1,
                                                   use_terminal_gap_penalty=False,
                                                   emulate_rb=False)
        assert sa == "gattacagattaca", "seqa incorrect"
        assert sb == "gattac---ttaca", "seqb incorrect"
        assert score == 76, "expected score=76"
        assert sts == "ok", "alignment failed"

    def test_align_it_aa_rb(self) -> None:
        """align_it_aa_rb() with test sequences, legacy (ruby) score matrix
        raises a ValidationError"""
        aligner = al.AlignIt()
        [sa, sb, score, sts] = aligner.align_it_aa(seqa="gattacagattaca",
                                                   seqb="gattacttaca", 
                                                   gap_ini=6,
                                                   gap_ext=1,
                                                   use_terminal_gap_penalty=False,
                                                   emulate_rb=True)
        assert sa == "gattacagattaca", "seqa incorrect"
        assert sb == "gattacttaca---", "seqb incorrect"
        assert score == 0, "expected score=0"
        assert sts == "ok", "alignment failed"