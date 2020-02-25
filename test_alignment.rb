require 'minitest/autorun'
require './alignment'

HXB2_NUCS = 'TGGAAGGGCTAATTCACTCCCAACGAAGACAAGATATCCTTGATCTGTGGATCTACCACACACAA'
DEFAULT_PENALTIES = [
    15, # gap open penalty
    5   # gap_extend_penalty
]

class TestAlignment < Minitest::Test
  def test_align_it_perfect_match
    hxb2 = HXB2_NUCS
    ref = seq = hxb2
    aligned_ref, aligned_seq = align_it(ref, seq, *DEFAULT_PENALTIES)
    assert_equal ref, aligned_ref
    assert_equal seq, aligned_seq
  end

  def test_align_it_mismatch
    ref = 'TGGAAGGGCTAATTCACTCCCAACGAAGACAAGATATCCTTGATCTGTGGATCTACCACACACAA'
    seq = 'TGGAAGGGATAATTCACTCCCAACGAAGACAAGATATCCTTGATCTGTGGATCTACCACACACAA'
    # mismatch:    ^
    aligned_ref, aligned_seq = align_it(ref, seq, *DEFAULT_PENALTIES)
    assert_equal ref, aligned_ref
    assert_equal seq, aligned_seq
  end

  def test_align_it_insert
    ref = \
        'TGGAAGGGCTAATTCACTCCCAACGAAGACAAGATATCCTTGATCTGTGGATCTACCACACA'
    seq = \
        'TGGAAGGGCTAATTCACTGAGACTCCCAACGAAGACAAGATATCCTTGATCTGTGGATCTACCACACA'
    # insertion:           ^^^^^^
    expected_aref = \
        '------TGGAAGGGCTAATTCACTCCCAACGAAGACAAGATATCCTTGATCTGTGGATCTACCACACA'
    gap_open_penalty = 15
    gap_extend_penalty = 5

    aligned_ref, aligned_seq = align_it(ref, seq, gap_open_penalty, gap_extend_penalty)
    assert_equal expected_aref, aligned_ref
    assert_equal seq, aligned_seq
  end

  def test_align_it_deletion
    ref = \
        'TGGAAGGGCTAATTCACTCCCAACGAAGACAAGATATCCTTGATCTGTGGATCTACCACACA'
    # deletion:            ^^^^^^
    seq = \
        'TGGAAGGGCTAATTCACTGAAGACAAGATATCCTTGATCTGTGGATCTACCACACA'
    expected_aseq = \
        '------TGGAAGGGCTAATTCACTGAAGACAAGATATCCTTGATCTGTGGATCTACCACACA'
    gap_open_penalty = 15
    gap_extend_penalty = 5

    aligned_ref, aligned_seq = align_it(ref, seq, gap_open_penalty, gap_extend_penalty)
    assert_equal ref, aligned_ref
    assert_equal expected_aseq, aligned_seq
  end

  def test_align_it_aa_match
    ref1 = 'R'
    seq1 = 'R'
    ref2 = 'D'
    seq2 = 'D'

    aref1, aseq1 = align_it_aa(ref1, seq1, *DEFAULT_PENALTIES)
    aref2, aseq2 = align_it_aa(ref2, seq2, *DEFAULT_PENALTIES)

    assert_equal [ref1, seq1], [aref1, aseq1]
    assert_equal [ref2, seq2], [aref2, aseq2]
  end

  def test_align_it_aa_mismatch
    ref = 'RRRRRRR'
    seq = 'RRRDRRR'
    # Main difference is that align_it_aa_rb uses same score for all combos.

    aref, aseq = align_it_aa(ref, seq, *DEFAULT_PENALTIES)

    assert_equal ref, aref
    assert_equal seq, aseq
  end

  #########################################
  # The following test the python versions of alignment ('gotoh').
  # alignment.cpp must be altered with the following changes:
  # - init_pairscore parameters must be changed
  # - for align_it_aa, init_pairscore must be replaced with init_pairscore_hiv25
  # - score must be returned
  # - degap calls must be removed

  # def test_align_it_perfect_match
  #   hxb2 = HXB2_NUCS
  #   ref = seq = hxb2
  #   aligned_ref, aligned_seq = align_it_py(ref, seq, *DEFAULT_PENALTIES)
  #   assert_equal aligned_ref, ref
  #   assert_equal aligned_seq, seq
  # end
  #
  # def test_align_it_mismatch
  #   ref = 'TGGAAGGGCTAATTCACTCCCAACGAAGACAAGATATCCTTGATCTGTGGATCTACCACACACAA'
  #   seq = 'TGGAAGGGATAATTCACTCCCAACGAAGACAAGATATCCTTGATCTGTGGATCTACCACACACAA'
  #   # mismatch:    ^
  #   expected_matches = ref.size - 1
  #   expected_mismatches = 1
  #   aligned_ref, aligned_seq, score = align_it_py(ref, seq, *DEFAULT_PENALTIES)
  #   assert_equal aligned_ref, ref
  #   assert_equal aligned_seq, seq
  #   assert_equal score, 5 * expected_matches - 4 * expected_mismatches
  # end
  #
  # def test_align_it_insert
  #   ref = \
  #       'TGGAAGGGCTAATTCACTCCCAACGAAGACAAGATATCCTTGATCTGTGGATCTACCACACA'
  #   seq = \
  #       'TGGAAGGGCTAATTCACTGAGACTCCCAACGAAGACAAGATATCCTTGATCTGTGGATCTACCACACA'
  #   # insertion:           ^^^^^^
  #   expected_aref = \
  #       'TGGAAGGGCTAATTCACT------CCCAACGAAGACAAGATATCCTTGATCTGTGGATCTACCACACA'
  #   expected_matches = ref.size
  #   expected_gap_count = 1
  #   expected_gap_size = 6
  #   gap_open_penalty = 15
  #   gap_extend_penalty = 5
  #   use_terminal_gap_penalty = 1
  #
  #   aligned_ref, aligned_seq, score = align_it_py(ref,
  #                                                 seq,
  #                                                 gap_open_penalty,
  #                                                 gap_extend_penalty,
  #                                                 use_terminal_gap_penalty)
  #   assert_equal aligned_ref, expected_aref
  #   assert_equal aligned_seq, seq
  #   assert_equal score, (5 * expected_matches -
  #       expected_gap_count * gap_open_penalty -
  #       expected_gap_size * gap_extend_penalty)
  # end
  #
  # def test_align_it_deletion
  #   ref = \
  #       'TGGAAGGGCTAATTCACTCCCAACGAAGACAAGATATCCTTGATCTGTGGATCTACCACACA'
  #   # deletion:            ^^^^^^
  #   seq = \
  #       'TGGAAGGGCTAATTCACTGAAGACAAGATATCCTTGATCTGTGGATCTACCACACA'
  #   expected_aseq = \
  #       'TGGAAGGGCTAATTCACT------GAAGACAAGATATCCTTGATCTGTGGATCTACCACACA'
  #   expected_matches = seq.size
  #   expected_gap_count = 1
  #   expected_gap_size = 6
  #   gap_open_penalty = 15
  #   gap_extend_penalty = 5
  #   use_terminal_gap_penalty = 1
  #
  #   aligned_ref, aligned_seq, score = align_it_py(ref,
  #                                                 seq,
  #                                                 gap_open_penalty,
  #                                                 gap_extend_penalty,
  #                                                 use_terminal_gap_penalty)
  #   assert_equal aligned_ref, ref
  #   assert_equal aligned_seq, expected_aseq
  #   assert_equal score, (5 * expected_matches -
  #       expected_gap_count * gap_open_penalty -
  #       expected_gap_size * gap_extend_penalty)
  # end
  #
  # def test_align_it_aa_match
  #   ref1 = 'R'
  #   seq1 = 'R'
  #   ref2 = 'D'
  #   seq2 = 'D'
  #   # align_it_aa_py uses different scores for each amino acid combination.
  #   expected_score1 = 7
  #   expected_score2 = 8
  #
  #   aref1, aseq1, score1 = align_it_aa_py(ref1, seq1, *DEFAULT_PENALTIES)
  #   aref2, aseq2, score2 = align_it_aa_py(ref2, seq2, *DEFAULT_PENALTIES)
  #
  #   assert_equal [aref1, aseq1], [ref1, seq1]
  #   assert_equal [aref2, aseq2], [ref2, seq2]
  #   assert_equal score1, expected_score1
  #   assert_equal score2, expected_score2
  # end
  #
  # def test_align_it_aa_py_mismatch
  #   ref0 = 'RRRRRRR'
  #   seq1 = 'RRRNRRR'
  #   seq2 = 'RRRDRRR'
  #   # align_it_aa_py uses different scores for each amino acid combination.
  #   expected_score1 = 6*7 - 5
  #   expected_score2 = 6*7 - 11
  #
  #   aref1, aseq1, score1 = align_it_aa_py(ref0, seq1, *DEFAULT_PENALTIES)
  #   aref2, aseq2, score2 = align_it_aa_py(ref0, seq2, *DEFAULT_PENALTIES)
  #
  #   assert_equal [aref1, aseq1], [ref0, seq1]
  #   assert_equal [aref2, aseq2], [ref0, seq2]
  #   assert_equal score1, expected_score1
  #   assert_equal score2, expected_score2
  # end

end