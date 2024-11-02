require_relative 'test_helper'
require_relative '../lib/cfe_gotoh'


class CfeGotohTest < Minitest::Test
end


class ScoreAlignmentTest < CfeGotohTest
  REGULAR_BASES = ['A', 'C', 'G', 'T'].freeze
  MIXTURES = {
    'R': ['A', 'G'],
    'Y': ['C', 'T'],
    'S': ['G', 'C'],
    'W': ['A', 'T'],
    'K': ['G', 'T'],
    'M': ['A', 'C'],
    'B': ['C', 'G', 'T'],
    'D': ['A', 'G', 'T'],
    'H': ['A', 'C', 'T'],
    'V': ['A', 'C', 'G'],
    'N': ['A', 'C', 'G', 'T']
  }.freeze

  def score_alignment_symmetric_test(expected, base1, base2)
    assert_equal expected, CfeGotoh.score_alignment(base1, base2)
    assert_equal expected, CfeGotoh.score_alignment(base2, base1)
  end

  def test_regular_base_scores
    REGULAR_BASES.each do |std_base|
      REGULAR_BASES.each do |query_base|
        expected = -1.0
        if (std_base == query_base)
          expected = 1.0
        end
        assert_equal expected, CfeGotoh.score_alignment(std_base, query_base)
      end
    end
  end

  def test_mixture_scores
    MIXTURES.each do |std_base, _|
      MIXTURES.each do |query_base, __|
        expected = -1.0
        if (std_base == query_base)
          expected = 1.0
        end
        assert_equal expected, CfeGotoh.score_alignment(std_base, query_base)
      end
    end
  end

  def test_mixture_with_regular_scores
    REGULAR_BASES.each do |base1|
      MIXTURES.each do |base2, possible_bases|
        expected = -1.0
        if (base2 == 'N')
          expected = -3.0
        elsif (possible_bases.include?(base1))
          expected = 1.0
        end
        score_alignment_symmetric_test(expected, base1, base2)
      end
    end 
  end

  def test_score_x
    (REGULAR_BASES + MIXTURES).each do |base|
      expected = -6.0
      if (base == 'N')
        expected = -1.0
      end
      score_alignment_symmetric_test(expected, base, 'X')
    end
  end

  def test_special_characters
    assert_equal 50.0, CfeGotoh.score_alignment('$', '$')
    score_alignment_symmetric_test(50.0, '$', '$')
    score_alignment_symmetric_test(1.0, 'T', 'U')
    assert_equal 0.0, CfeGotoh.score_alignment('N', 'N')
    score_alignment_symmetric_test(3.0, 'X', '-')
    REGULAR_BASES.each do |base|
      score_alignment_symmetric_test(1.0, base, '*')
      score_alignment_symmetric_test(0.7, base, '&')
      score_alignment_symmetric_test(0.0, base, '$')
      score_alignment_symmetric_test(-20.0, base, '.')
    end
  end

  def test_multiple_bases
    seq1 = 'AGACTCTVC---'
    seq2 = 'CGANNCTGCXXX'
    expected = -1.0 + 2.0 - 3.0 - 3.0 + 4.0 + 9.0
    assert_equal expected, CfeGotoh.score_alignment(seq1, seq2)
  end

  MAKE_GAP_LIST_TEST_CASES = [
    {
      name: 'empty_sequence',
      seq: '',
      expected: []
    },
    {
      name: 'no_gaps',
      seq: 'ACAGAT',
      expected: []
    },
    {
      name: 'gap_in_middle',
      seq: 'ACA--GATC',
      expected: [[3, 4]]
    },
    {
      name: 'gap_at_start',
      seq: '--ACAGATCC',
      expected: [[0, 1]]
    },
    {
      name: 'gap_at_end',
      seq: 'ACAGATCC----',
      expected: [[8, 9, 10, 11]]
    },
    {
      name: 'multiple_gaps_in_middle',
      seq: 'ACA--GA-TC',
      expected: [[3, 4], [7]]
    },
    {
      name: 'multiple_gaps_throughout',
      seq: '-ACA--GA-TC----',
      expected: [[0], [4, 5], [8], [11, 12, 13, 14]]
    }
  ]
end
