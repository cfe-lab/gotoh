require 'minitest/autorun'

require './alignment'

HXB2_NUCS = 'TGGAAGGGCTAATTCACTCCCAACGAAGACAAGATATCCTTGATCTGTGGATCTACCACACACAA'
DEFAULT_PENALTIES = [
    15,  # gap open penalty
    5,  # gap extend penalty
]
class TestAlginment < Minitest::Test
  def test_align_it_perfect_match
    hxb2 = HXB2_NUCS
    ref = seq = hxb2
    aligned_ref, aligned_seq = align_it(ref, seq, *DEFAULT_PENALTIES)
    assert_equal(aligned_ref, ref)
    assert_equal(aligned_seq, seq)
  end
end