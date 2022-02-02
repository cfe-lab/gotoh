"""Test the align_it module"""

import alignment.align_it as al
import pytest
from pydantic.error_wrappers import ValidationError


HXB2_NUCS = 'TGGAAGGGCTAATTCACTCCCAACGAAGACAAGATATCCTTGATCTGTGGATCTACCACACACAA'
GIP = 15
GEP = 5
USE_TERM_GAP = True

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


    def test_align_it_perfect_match(self):
        hxb2 = HXB2_NUCS
        ref = seq = hxb2
        aligner = al.AlignIt()
        [aref, aseq, score, sts] = aligner.align_it(ref, seq,
                                                    GIP, GEP, USE_TERM_GAP)
        assert aref == hxb2
        assert aseq == hxb2
        assert score == 5 * len(hxb2)
        assert sts == "ok", "alignment failed"


    def test_align_it_mismatch(self):
        ref = 'TGGAAGGGCTAATTCACTCCCAACGAAGACAAGATATCCTTGATCTGTGGATCTACCACACACAA'
        seq = 'TGGAAGGGATAATTCACTCCCAACGAAGACAAGATATCCTTGATCTGTGGATCTACCACACACAA'
        # mismatch:    ^
        aligner = al.AlignIt()
        expected_matches = len(ref) - 1
        expected_mismatches = 1
        [aref, aseq, score, sts] = aligner.align_it(ref, seq,
                                                    GIP, GEP, USE_TERM_GAP)
        assert aref == ref
        assert aseq == seq
        assert score == 5 * expected_matches - 4 * expected_mismatches
        assert sts == "ok", "alignment failed"


    # noinspection DuplicatedCode
    def test_align_it_insert(self):
        aligner = al.AlignIt()
        ref = \
            'TGGAAGGGCTAATTCACTCCCAACGAAGACAAGATATCCTTGATCTGTGGATCTACCACACA'
        seq = \
            'TGGAAGGGCTAATTCACTGAGACTCCCAACGAAGACAAGATATCCTTGATCTGTGGATCTACCACACA'
        # insertion:           ^^^^^^
        expected_aref = \
            'TGGAAGGGCTAATTCACT------CCCAACGAAGACAAGATATCCTTGATCTGTGGATCTACCACACA'
        expected_matches = len(ref)
        expected_gap_count = 1
        expected_gap_size = 6
        gap_open_penalty = 15
        gap_extend_penalty = 5
        use_terminal_gap_penalty = 1

        [aref, aseq, score, sts] = aligner.align_it(ref,
                                                    seq,
                                                    gap_open_penalty,
                                                    gap_extend_penalty,
                                                    use_terminal_gap_penalty)
        assert aref == expected_aref
        assert aseq == seq
        assert score == (5 * expected_matches -
                        expected_gap_count * gap_open_penalty -
                        expected_gap_size * gap_extend_penalty)
        assert sts == "ok", "alignment failed"


    # noinspection DuplicatedCode
    def test_align_it_deletion(self):
        aligner = al.AlignIt()
        ref = \
            'TGGAAGGGCTAATTCACTCCCAACGAAGACAAGATATCCTTGATCTGTGGATCTACCACACA'
        # deletion:            ^^^^^^
        seq = \
            'TGGAAGGGCTAATTCACTGAAGACAAGATATCCTTGATCTGTGGATCTACCACACA'
        expected_aseq = \
            'TGGAAGGGCTAATTCACT------GAAGACAAGATATCCTTGATCTGTGGATCTACCACACA'
        expected_matches = len(seq)
        expected_gap_count = 1
        expected_gap_size = 6
        gap_open_penalty = 15
        gap_extend_penalty = 5
        use_terminal_gap_penalty = 1

        [aref, aseq, score, sts] = aligner.align_it(ref,
                                                    seq,
                                                    gap_open_penalty,
                                                    gap_extend_penalty,
                                                    use_terminal_gap_penalty)
        assert aref == ref
        assert aseq == expected_aseq
        assert score == (5 * expected_matches -
                        expected_gap_count * gap_open_penalty -
                        expected_gap_size * gap_extend_penalty)
        assert sts == "ok", "alignment failed"


    # noinspection DuplicatedCode
    def test_align_it_aa_match(self):
        aligner = al.AlignIt()
        ref1 = 'R'
        seq1 = 'R'
        ref2 = 'D'
        seq2 = 'D'
        # align_it_aa uses different scores for each amino acid combination.
        expected_score1 = 7
        expected_score2 = 8

        [aref1, aseq1, score1, sts] = aligner.align_it_aa(ref1, seq1,
                                                          GIP, GEP, USE_TERM_GAP)
        [aref2, aseq2, score2, sts] = aligner.align_it_aa(ref2, seq2,
                                                           GIP, GEP, USE_TERM_GAP)

        assert (aref1, aseq1) == (ref1, seq1)
        assert (aref2, aseq2) == (ref2, seq2)
        assert score1 == expected_score1
        assert score2 == expected_score2
        assert sts == "ok", "alignment failed"


    # noinspection DuplicatedCode
    def test_align_it_aa_mismatch(self):
        aligner = al.AlignIt()
        ref0 = 'RRRRRRR'
        seq1 = 'RRRNRRR'
        seq2 = 'RRRDRRR'
        # align_it_aa uses different scores for each amino acid combination.
        expected_score1 = 6*7 - 5
        expected_score2 = 6*7 - 11

        [aref1, aseq1, score1, sts] = aligner.align_it_aa(ref0, seq1,
                                                          GIP, GEP, USE_TERM_GAP)
        [aref2, aseq2, score2, sts] = aligner.align_it_aa(ref0, seq2,
                                                          GIP, GEP, USE_TERM_GAP)

        assert (aref1, aseq1) == (ref0, seq1)
        assert (aref2, aseq2) == (ref0, seq2)
        assert score1 == expected_score1
        assert score2 == expected_score2
        assert sts == "ok", "alignment failed"


    # noinspection DuplicatedCode
    def test_align_it_aa_rb_mismatch(self):
        aligner = al.AlignIt()
        ref = 'RRRRRRR'
        seq = 'RRRDRRR'
        # Main difference is that align_it_aa_rb uses same score for all combos.

        [aref, aseq, _, sts] = aligner.align_it_aa(ref, seq, 
                                                   GIP, GEP, 
                                                   emulate_rb=True)

        assert aref == ref
        assert aseq == seq
        assert sts == "ok", "alignment failed"
