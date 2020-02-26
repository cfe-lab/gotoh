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
        'TGGAAGGGCTAATTCACT------CCCAACGAAGACAAGATATCCTTGATCTGTGGATCTACCACACA'
    gap_open_penalty = 3
    gap_extend_penalty = 1

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
        'TGGAAGGGCTAATTCACT------GAAGACAAGATATCCTTGATCTGTGGATCTACCACACA'
    gap_open_penalty = 3
    gap_extend_penalty = 1

    aligned_ref, aligned_seq = align_it(ref, seq, gap_open_penalty, gap_extend_penalty)
    assert_equal ref, aligned_ref
    assert_equal expected_aseq, aligned_seq
  end

  def test_align_it_aa_mismatch
    ref = 'RRRRRRR'
    seq = 'RRRDRRR'

    aref, aseq = align_it_aa(ref, seq, *DEFAULT_PENALTIES)

    assert_equal ref, aref
    assert_equal seq, aseq
  end

end