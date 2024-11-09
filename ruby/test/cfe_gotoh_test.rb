require_relative 'test_helper'
require_relative '../lib/cfe_gotoh'


class CfeGotohTest < Minitest::Test
end


class ScoreAlignmentTest < CfeGotohTest
  REGULAR_BASES = ['A', 'C', 'G', 'T'].freeze
  MIXTURES = {
    'R' => ['A', 'G'],
    'Y' => ['C', 'T'],
    'S' => ['G', 'C'],
    'W' => ['A', 'T'],
    'K' => ['G', 'T'],
    'M' => ['A', 'C'],
    'B' => ['C', 'G', 'T'],
    'D' => ['A', 'G', 'T'],
    'H' => ['A', 'C', 'T'],
    'V' => ['A', 'C', 'G'],
    'N' => ['A', 'C', 'G', 'T']
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
        score_alignment_symmetric_test(expected, std_base, query_base)
      end
    end
  end

  def test_mixture_scores
    MIXTURES.each do |std_base, _|
      MIXTURES.each do |query_base, __|
        expected = -1.0
        if (std_base == query_base)
          expected = 1.0
          if (std_base == 'N')
            expected = 0.0
          end
        end
        score_alignment_symmetric_test(expected, std_base, query_base)
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
    (REGULAR_BASES + MIXTURES.keys()).each do |base|
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
end


class MakeGapListTest < CfeGotohTest
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

  MAKE_GAP_LIST_TEST_CASES.each do |test_entry|
    define_method("test_#{test_entry[:name]}") do
      assert_equal test_entry[:expected], CfeGotoh.make_gap_list(test_entry[:seq])
    end
  end
end


class TrimLeadingDashesTest < CfeGotohTest
  TRIM_LEADING_DASHES_TEST_CASES = [
    {
      name: 'no_leading_dashes',
      std: 'ACAGAT',
      query: 'ACACAT',
      expected_std: 'ACAGAT',
      expected_query: 'ACACAT'
    },
    {
      name: 'one_leading_dash',
      std: '-ACAGAT',
      query: 'GACACAT',
      expected_std: 'ACAGAT',
      expected_query: 'ACACAT'
    },
    {
      name: 'several_leading_dashes',
      std: '----ACAGAT',
      query: 'GGGGACACAT',
      expected_std: 'ACAGAT',
      expected_query: 'ACACAT'
    },
    {
      name: 'no_leading_dashes_other_dashes_ignored',
      std: 'ACA---CATGAT-',
      query: 'ACAGGG---CATC',
      expected_std: 'ACA---CATGAT-',
      expected_query: 'ACAGGG---CATC'
    },
    {
      name: 'one_leading_dash_other_dashes_ignored',
      std: '-ACA---CATGAT-',
      query: '-ACAGGG---CATC',
      expected_std: 'ACA---CATGAT-',
      expected_query: 'ACAGGG---CATC'
    },
    {
      name: 'several_leading_dashes_other_dashes_ignored',
      std: '----ACA---CATGAT-',
      query: 'GGGGACAGGG---CATC',
      expected_std: 'ACA---CATGAT-',
      expected_query: 'ACAGGG---CATC'
    }
  ]

  TRIM_LEADING_DASHES_TEST_CASES.each do |test_entry|
    define_method("test_#{test_entry[:name]}") do
      std = test_entry[:std]
      query = test_entry[:query]
      CfeGotoh.trim_leading_dashes(std, query)
      assert_equal test_entry[:expected_std], std
      assert_equal test_entry[:expected_query], query
    end
  end
end


class TrimTrailingDashesTest < CfeGotohTest
  TRIM_TRAILING_DASHES_TEST_CASES = [
    {
      name: 'no_trailing_dashes',
      std: 'ACAGAT',
      query: 'ACACAT',
      expected_std: 'ACAGAT',
      expected_query: 'ACACAT'
    },
    {
      name: 'one_trailing_dash',
      std: 'ACAGAT-',
      query: 'ACACATG',
      expected_std: 'ACAGAT',
      expected_query: 'ACACAT'
    },
    {
      name: 'several_trailing_dashes',
      std: 'ACAGAT----',
      query: 'ACACATGGGG',
      expected_std: 'ACAGAT',
      expected_query: 'ACACAT'
    },
    {
      name: 'no_trailing_dashes_other_dashes_ignored',
      std: '-ACA---CATGAT',
      query: 'CACAGGG---CAT',
      expected_std: '-ACA---CATGAT',
      expected_query: 'CACAGGG---CAT'
    },
    {
      name: 'one_trailing_dash_other_dashes_ignored',
      std: '-ACA---CATGAT-',
      query: 'CACAGGG---CATC',
      expected_std: '-ACA---CATGAT',
      expected_query: 'CACAGGG---CAT'
    },
    {
      name: 'several_trailing_dashes_other_dashes_ignored',
      std: '-ACA---CATGAT----',
      query: 'CACAGGG---CATCGGC',
      expected_std: '-ACA---CATGAT',
      expected_query: 'CACAGGG---CAT'
    }
  ]

  TRIM_TRAILING_DASHES_TEST_CASES.each do |test_entry|
    define_method("test_#{test_entry[:name]}") do
      std = test_entry[:std]
      query = test_entry[:query]
      CfeGotoh.trim_trailing_dashes(std, query)
      assert_equal test_entry[:expected_std], std
      assert_equal test_entry[:expected_query], query
    end
  end
end


class FixIncompleteEdgeCodonTest < CfeGotohTest
  FIX_INCOMPLETE_EDGE_CODON_TEST_CASES = [
    {
      name: 'no_leading_dashes',
      seq: 'ACTAGG',
      expected: 'ACTAGG'
    },
    {
      name: 'one_complete_blank_codon_leading',
      seq: '---ACTAGG',
      expected: '---ACTAGG'
    },
    {
      name: 'several_complete_blank_codons_leading',
      seq: '---------ACTAGG',
      expected: '---------ACTAGG'
    },
    {
      name: 'no_leading_dashes_other_blanks_ignored',
      seq: 'ACT---AGG---',
      expected: 'ACT---AGG---'
    },
    {
      name: 'one_leading_dash',
      seq: '-GGACTAGG',
      expected: '---ACTAGG'
    },
    {
      name: 'two_leading_dashes',
      seq: '--GACTAGG',
      expected: '---ACTAGG'
    },
    {
      name: 'one_dash_plus_full_codon_leading',
      seq: '----GGACTAGG',
      expected: '------ACTAGG'
    },
    {
      name: 'two_dashes_plus_full_codons_leading',
      seq: '--------GACTAGG',
      expected: '---------ACTAGG'
    },
    {
      name: 'leading_dashes_other_blanks_ignored',
      seq: '-----GACT---AGG---',
      expected: '------ACT---AGG---'
    },
    {
      name: 'no_trailing_dashes',
      seq: 'ACTAGG',
      expected: 'ACTAGG',
      side: :trailing
    },
    {
      name: 'one_complete_blank_codon_trailing',
      seq: 'ACTAGG---',
      expected: 'ACTAGG---',
      side: :trailing
    },
    {
      name: 'several_complete_blank_codons_trailing',
      seq: 'ACTAGG---------',
      expected: 'ACTAGG---------',
      side: :trailing
    },
    {
      name: 'no_trailing_dashes_other_blanks_ignored',
      seq: '---ACT---AGG',
      expected: '---ACT---AGG',
      side: :trailing
    },
    {
      name: 'one_trailing_dash',
      seq: 'ACTAGGTT-',
      expected: 'ACTAGG---',
      side: :trailing
    },
    {
      name: 'two_trailing_dashes',
      seq: 'ACTAGGT--',
      expected: 'ACTAGG---',
      side: :trailing
    },
    {
      name: 'one_trailing_dash_plus_full_codon',
      seq: 'ACTAGGTT----',
      expected: 'ACTAGG------',
      side: :trailing
    },
    {
      name: 'two_trailing_dashes_plus_full_codons',
      seq: 'ACTAGGT--------',
      expected: 'ACTAGG---------',
      side: :trailing
    },
    {
      name: 'trailing_dashes_other_blanks_ignored',
      seq: '------ACT---AGGT-----',
      expected: '------ACT---AGG------',
      side: :trailing
    }
  ]

  FIX_INCOMPLETE_EDGE_CODON_TEST_CASES.each do |test_entry|
    define_method("test_#{test_entry[:name]}") do
      seq = test_entry[:seq]
      side = test_entry[:side]
      if (side.nil?)
        side = :leading
      end
      CfeGotoh.fix_incomplete_edge_codon(seq, side)
      assert_equal test_entry[:expected], seq
    end
  end
end


# class merge_insertions_and_deletions_to_fix_oof_sequences
class MergeInsertionsAndDeletionsToFixOofSequencesTest < CfeGotohTest
  def test_standard_and_query_must_be_same_length
    assert_raises RuntimeError do
      CfeGotoh.merge_insertions_and_deletions_to_fix_oof_sequences('ACT', 'ACTACT')
    end
  end

  MERGE_INDELS_TEST_CASES = [
    {
      name: 'no_need_for_merge',
      std: 'ACTAAG',
      query: 'ACTCAG',
      expected_std: 'ACTAAG',
      expected_query: 'ACTCAG'
    },
    {
      name: 'blanks_but_no_need_for_merge',
      std: 'ACT---AGCTTT',
      query: 'ACTAGC---TTC',
      expected_std: 'ACT---AGCTTT',
      expected_query: 'ACTAGC---TTC'
    },
    {
      name: 'bad_length_but_no_possible_merges',
      std: 'ACTAAGC',
      query: 'ACTCAGC',
      expected_std: 'ACTAAGC',
      expected_query: 'ACTCAGC'
    },
    {
      name: 'merge_one_base_ahead',
      std: 'ACT-AAG',
      query: 'ACTC-AG',
      expected_std: 'ACTAAG',
      expected_query: 'ACTCAG'
    },
    {
      name: 'merge_one_base_behind',
      std: 'ACT-AAG',
      query: 'AC-TCAG',
      expected_std: 'ACTAAG',
      expected_query: 'ACTCAG'
    },
    {
      name: 'merge_two_bases_ahead',
      std: 'ACT-AAG',
      query: 'ACTCA-G',
      expected_std: 'ACTAAG',
      expected_query: 'ACTCAG'
    },
    {
      name: 'merge_two_bases_behind',
      std: 'ACT-AAG',
      query: 'A-CTCAG',
      expected_std: 'ACTAAG',
      expected_query: 'ACTCAG'
    },
    {
      name: 'one_base_behind_preferred',
      std: 'ACT-AAG',
      query: 'AC-T-AG',
      expected_std: 'ACTAAG',
      expected_query: 'ACT-AG'
    },
    {
      name: 'one_base_ahead_preferred_over_two_bases_behind',
      std: 'ACT-AAA',
      query: 'A-CT-AA',
      expected_std: 'ACTAAA',
      expected_query: 'A-CTAA'
    },
    {
      name: 'two_bases_behind_preferred_over_two_bases_ahead',
      std: 'ACT-AAA',
      query: 'A-CTA-A',
      expected_std: 'ACTAAA',
      expected_query: 'ACTA-A'
    },
    {
      name: 'gaps_too_far_to_merge',
      std: 'ACT-AAA',
      query: 'ACTAAA-',
      expected_std: 'ACT-AAA',
      expected_query: 'ACTAAA-'
    },
    {
      name: 'merges_stop_at_cogent_length',
      std: 'ACT-AAAG-GG-CC',
      query: 'AC-TAAAGGG--CC',
      expected_std: 'ACTAAAGGG-CC',
      expected_query: 'ACTAAAGGG-CC'
    },
    {
      name: 'impossible_merge_skipped_but_later_ones_happen',
      std: 'ACT-AAAGGG-CC',
      query: 'ACTAAA-GG-GCC',
      expected_std: 'ACT-AAAGGGCC',
      expected_query: 'ACTAAA-GGGCC'
    }
  ]

  MERGE_INDELS_TEST_CASES.each do |test_entry|
    define_method("test_#{test_entry[:name]}") do
      std = test_entry[:std]
      query = test_entry[:query]
      side = test_entry[:side]
      CfeGotoh.merge_insertions_and_deletions_to_fix_oof_sequences(std, query)
      assert_equal test_entry[:expected_std], std
      assert_equal test_entry[:expected_query], query
    end
  end
end


class ClusterGapsTest < CfeGotohTest
  CLUSTER_GAPS_TEST_CASES = [
    {
      name: 'no_gaps',
      gaps: [],
      expected: []
    },
    {
      name: 'empty_gap_causes_nothing',
      gaps: [[]],
      expected: []
    },
    {
      name: 'good_size_gap',
      gaps: [[3, 4, 5]],
      expected: [[3, 4, 5]]
    },
    {
      name: 'bad_size_gap',
      gaps: [[3, 4, 5, 6]],
      expected: [[3, 4, 5, 6]]
    },
    {
      name: 'merge_two_close_gaps_to_first',
      gaps: [[10, 11], [13]],
      expected: [[10, 11, 12]]
    },
    {
      name: 'merge_two_close_gaps_to_second',
      gaps: [[7], [10, 11]],
      expected: [[9, 10, 11]]
    },
    {
      name: 'two_gaps_can_merge_edge_case',
      gaps: [[2, 3, 4, 5, 6], [14]],
      expected: [[2, 3, 4, 5, 6, 7]]
    },
    {
      name: 'two_gaps_too_far_edge_case',
      gaps: [[2, 3, 4, 5, 6], [15]],
      expected: [[2, 3, 4, 5, 6], [15]]
    },
    {
      name: 'two_gaps_too_far',
      gaps: [[2, 3, 4, 5, 6], [21]],
      expected: [[2, 3, 4, 5, 6], [21]]
    },
    {
      name: 'three_close_gaps_merge',
      gaps: [[8], [12, 13, 14, 15], [18, 19, 20, 21]],
      expected: [[11, 12, 13, 14, 15, 16, 17, 18, 19]]
    },
    {
      name: 'three_gaps_merge_edge_case',
      gaps: [[8], [12, 13, 14, 15], [19, 20, 21, 22]],
      expected: [[11, 12, 13, 14, 15, 16, 17, 18, 19]]
    },
    {
      name: 'three_gaps_too_far_edge_case',
      gaps: [[8], [12, 13, 14, 15], [20, 21, 22, 23]],
      expected: [[8], [12, 13, 14, 15], [20, 21, 22, 23]]
    },
    {
      name: 'typical_case',
      gaps: [[3, 4, 5], [8, 9], [13], [19, 20, 21, 22, 23, 24], [27], [32], [38, 39, 40, 41], [50, 51], [60], [70], [75, 76]],
      expected: [[3, 4, 5], [8, 9, 10], [19, 20, 21, 22, 23, 24], [31, 32, 33, 34, 35, 36], [50, 51], [60], [74, 75, 76]]
    }
  ]

  CLUSTER_GAPS_TEST_CASES.each do |test_entry|
    define_method("test_#{test_entry[:name]}") do
      assert_equal test_entry[:expected], CfeGotoh.cluster_gaps(test_entry[:gaps])
    end
  end

  def test_bad_gap_causes_error
    assert_raises CfeGotoh::GapMergeError do
      CfeGotoh.cluster_gaps([[3, 4]], raise_errors=true)
    end
  end

  def test_bad_gap_among_several_gaps_causes_error
    assert_raises CfeGotoh::GapMergeError do
      CfeGotoh.cluster_gaps([[3, 4, 5], [9, 10, 11, 12, 13, 14], [17]], raise_errors=true)
    end
  end
end


class AlignGapsToFrameTest < CfeGotohTest
  NO_COMMON_POSITIONS_TEST_CASES = [
    {
      name: 'no_gaps',
      gaps: [],
      expected: []
    },
    {
      name: 'already_in_frame',
      gaps: [[6, 7, 8]],
      expected: [[6, 7, 8]]
    },
    {
      name: 'shift_toward_beginning',
      gaps: [[7, 8, 9]],
      expected: [[6, 7, 8]]
    },
    {
      name: 'shift_toward_end',
      gaps: [[5, 6, 7]],
      expected: [[6, 7, 8]]
    },
    {
      name: 'two_in_frame',
      gaps: [[6, 7, 8], [15, 16, 17]],
      expected: [[6, 7, 8], [15, 16, 17]]
    },
    {
      name: 'two_needing_shifts',
      gaps: [[5, 6, 7, 8, 9, 10], [16, 17, 18]],
      expected: [[6, 7, 8, 9, 10, 11], [15, 16, 17]]
    },
    {
      name: 'several_in_frame',
      gaps: [[6, 7, 8], [15, 16, 17], [24, 25, 26]],
      expected: [[6, 7, 8], [15, 16, 17], [24, 25, 26]]
    },
    {
      name: 'several_needing_shifts',
      gaps: [[5, 6, 7, 8, 9, 10], [16, 17, 18], [31, 32, 33]],
      expected: [[6, 7, 8, 9, 10, 11], [15, 16, 17], [30, 31, 32]]
    }
  ]

  NO_COMMON_POSITIONS_TEST_CASES.each do |test_entry|
    define_method("test_#{test_entry[:name]}") do
      assert_equal test_entry[:expected], CfeGotoh.align_gaps_to_frame(test_entry[:gaps])
    end
  end

  WITH_COMMON_POSITIONS_TEST_CASES = [
    {
      name: 'no_gaps_with_common',
      gaps: [],
      common: [7, 15],
      expected: []
    },
    {
      name: 'too_far_from_common',
      gaps: [[3, 4, 5]],
      common: [7],
      expected: [[3, 4, 5]]
    },
    {
      name: 'too_far_from_common_but_needs_shift',
      gaps: [[2, 3, 4]],
      common: [7],
      expected: [[3, 4, 5]]
    },
    {
      name: 'too_far_before_common_edge_case',
      gaps: [[9, 10, 11]],
      common: [7],
      expected: [[9, 10, 11]]
    },
    {
      name: 'within_range_before_common_edge_case',
      gaps: [[12, 13, 14]],
      common: [7],
      expected: [[21, 22, 23]]
    },
    {
      name: 'within_range_after_common_edge_case',
      gaps: [[30, 31, 32]],
      common: [7],
      expected: [[21, 22, 23]]
    },
    {
      name: 'too_far_after_common_edge_case',
      gaps: [[33, 34, 35]],
      common: [7],
      expected: [[33, 34, 35]]
    },
    {
      name: 'too_far_after_common',
      gaps: [[45, 46, 47]],
      common: [7],
      expected: [[45, 46, 47]]
    },
    {
      name: 'too_far_after_common_but_needs_shift',
      gaps: [[46, 47, 48]],
      common: [7],
      expected: [[45, 46, 47]]
    },
    {
      name: 'offset_from_first_shifted_to_common_is_factored_into_second_shifted_to_common',
      # [57, 58, 59] is at codon 20 of the aligned sequence, which would be
      # at position 54 of the "raw" sequence; this should be just in the 
      # "catchment area" of the common insertion at codon 15.
      gaps: [[15, 16, 17], [57, 58, 59]],
      common: [7, 15],
      expected: [[21, 22, 23], [48, 49, 50]]
    },
    {
      name: 'offset_from_first_not_shifted_to_common_is_factored_into_second_shifted_to_common',
      gaps: [[4, 5, 6], [57, 58, 59]],
      common: [15],
      expected: [[3, 4, 5], [48, 49, 50]]
    },
    {
      name: 'offset_from_first_without_shifting_is_factored_into_second_shifted_to_common',
      gaps: [[3, 4, 5], [57, 58, 59]],
      common: [15],
      expected: [[3, 4, 5], [48, 49, 50]]
    },
    {
      name: 'offset_from_first_shifted_to_common_is_factored_into_second',
      # Even though [36, 37, 38] is in the "catchment area" of the common
      # insertion at codon 15, that's in the coordinates of the aligned
      # sequence; when the offset is accounted for, it should not be shifted.
      gaps: [[15, 16, 17], [36, 37, 38]],
      common: [7, 15],
      expected: [[21, 22, 23], [36, 37, 38]]
    },
    {
      name: 'offset_from_first_not_shifted_to_common_is_factored_into_second',
      gaps: [[3, 4, 5], [36, 37, 38]],
      common: [7, 15],
      expected: [[3, 4, 5], [36, 37, 38]]
    },
    {
      name: 'offset_from_first_already_at_common_is_factored_into_second',
      gaps: [[21, 22, 23], [36, 37, 38]],
      common: [7, 15],
      expected: [[21, 22, 23], [36, 37, 38]]
    },
    {
      name: 'offsets_taken_into_account_in_shifting_to_common',
      gaps: [[3, 4, 5, 6, 7, 8], [36, 37, 38], [111, 112, 113]],
      common: [15, 31],
      # [102, 103, 104] is codon 31 after the 9 base offset
      expected: [[3, 4, 5, 6, 7, 8], [36, 37, 38], [102, 103, 104]]
    },
    {
      name: 'two_gaps_shifted_to_same_common_position',
      gaps: [[14, 15, 16], [22, 23, 24]],
      common: [7],
      expected: [[21, 22, 23], [24, 25, 26]]
    },
    {
      name: 'typical_case',
      gaps: [[3, 4, 5], [17, 18, 19, 20, 21, 22], [40, 41, 42], [45, 46, 47]],
      common: [6, 14],
      expected: [[3, 4, 5], [21, 22, 23, 24, 25, 26], [39, 40, 41], [54, 55, 56]]
    }
  ]

  WITH_COMMON_POSITIONS_TEST_CASES.each do |test_entry|
    define_method("test_#{test_entry[:name]}") do
      assert_equal(
        test_entry[:expected],
        CfeGotoh.align_gaps_to_frame(test_entry[:gaps], test_entry[:common])
      )
    end
  end
end


class SpliceGapsIntoSequenceTest < CfeGotohTest
  SPLICE_GAPS_TEST_CASES = [
    {
      name: 'no_gaps',
      seq: 'ACTAAG',
      gaps: [],
      expected: 'ACTAAG'
    },
    {
      name: 'no_sequence',
      seq: '',
      gaps: [],
      expected: ''
    },
    {
      name: 'insert_at_beginning',
      seq: '---ACTAAG',
      gaps: [[0, 1, 2]],
      expected: '---ACTAAG'
    },
    {
      name: 'insert_at_beginning_removed_from_elsewhere',
      seq: 'ACT---AAG',
      gaps: [[0, 1, 2]],
      expected: '---ACTAAG'
    },
    {
      name: 'insert_at_end',
      seq: 'ACTAAG---',
      gaps: [[6, 7, 8]],
      expected: 'ACTAAG---'
    },
    {
      name: 'insert_at_end_removed_from_elsewhere',
      seq: 'ACT---AAG',
      gaps: [[6, 7, 8]],
      expected: 'ACTAAG---'
    },
    {
      name: 'insert_in_middle_retained',
      seq: 'ACT---AAG',
      gaps: [[3, 4, 5]],
      expected: 'ACT---AAG'
    },
    {
      name: 'insert_in_middle_corrected',
      seq: 'AC--TA-AG',
      gaps: [[3, 4, 5]],
      expected: 'ACT---AAG'
    },
    {
      name: 'offsets_accounted_for',
      seq: 'AACAT---GGG---G',
      gaps: [[3, 4, 5], [9, 10, 11]],
      expected: 'AAC---ATG---GGG'
    },
    {
      name: 'typical_case',
      seq: '---AACAT---GGG---G------',
      gaps: [[0, 1, 2], [6, 7, 8], [12, 13, 14], [18, 19, 20]],
      expected: '---AAC---ATG---GGG---'
    }
  ]

  SPLICE_GAPS_TEST_CASES.each do |test_entry|
    define_method("test_#{test_entry[:name]}") do
      assert_equal(
        test_entry[:expected],
        CfeGotoh.splice_gaps_into_sequence(test_entry[:seq], test_entry[:gaps])
      )
    end
  end
end


class FrameAlignTest < CfeGotohTest
  def test_bad_inserted_bases_error
    std = 'ACGTACGT-ACGT'
    query = 'ACGTACGTAACGT'
    assert_raises RuntimeError do
      CfeGotoh.frame_align(std, query, 3, 1, nil, false, true, true)
    end
  end

  def test_bad_inserted_bases_error
    std = 'ACGTACGTAACGT'
    query = 'ACGTACGT-ACGT'
    assert_raises RuntimeError do
      CfeGotoh.frame_align(std, query, 3, 1, nil, false, true, true)
    end
  end

  def test_unmerged_inserts_raise_error
    std = 'ACG--TACGTACGTAC-GT'
    query = 'ACGGGTACGTACGTACCGT'
    assert_raises RuntimeError do
      CfeGotoh.frame_align(std, query, 3, 1, nil, false, true, true)
    end
  end

  def test_unmerged_deletions_raise_error
    std = 'ACGGGTACGTACGTACCGT'
    query = 'ACG--TACGTACGTAC-GT'
    assert_raises RuntimeError do
      CfeGotoh.frame_align(std, query, 3, 1, nil, false, true, true)
    end
  end

  PREALIGN_TEST_CASES = [
    {
      name: 'edges_are_trimmed',
      std: '------ACGTACGTACGT------',
      query: 'AAAAAA-CGTACGTAC--AAAAAA',
      expected_std: 'ACGTACGTACGT',
      expected_query: '---TACGTA---'
    },
    {
      name: 'edges_not_trimmed_when_not_specified',
      std: '------ACGTACGTACGT------',
      query: 'AAAAAAACGTACGTACGTAAAAAA',
      trim: false,
      expected_std: '------ACGTACGTACGT------',
      expected_query: 'AAAAAAACGTACGTACGTAAAAAA'
    },
    {
      name: 'indels_are_merged',
      std: 'ACGT-ACGTACGT',
      query: 'ACGTAC-GTACGT',
      expected_std: 'ACGTACGTACGT',
      expected_query: 'ACGTACGTACGT'
    },
    {
      name: 'insertions_are_clustered',
      std: 'ACG--TA-CGTACGT',
      query: 'ACGGGTAACGTACGT',
      expected_std: 'ACG---TACGTACGT',
      expected_query: 'ACGGGTAACGTACGT'
    },
    {
      name: 'deletions_are_clustered',
      std: 'ACGGGTAACGTACGT',
      query: 'ACG--T-ACGTACGT',
      expected_std: 'ACGGGTAACGTACGT',
      expected_query: 'ACG---TACGTACGT'
    },
    {
      name: 'insertions_are_frame_aligned',
      std: 'ACGT---ACGTACGT',
      query: 'ACGTTTTACGTACGT',
      expected_std: 'ACG---TACGTACGT',
      expected_query: 'ACGTTTTACGTACGT'
    },
    {
      name: 'deletions_are_frame_aligned',
      std: 'ACGTTTTACGTACGT',
      query: 'ACGT---ACGTACGT',
      expected_std: 'ACGTTTTACGTACGT',
      expected_query: 'ACG---TACGTACGT'
    },
    {
      name: 'insertions_moved_to_common_positions',
      std: 'ACGTTTTACGTACGT',
      query: 'ACG---TACGTACGT',
      common_insert_locations: [3],
      expected_std: 'ACGTTTTACGTACGT',
      expected_query: 'ACGTACGTA---CGT',
    }
  ]

  PREALIGN_TEST_CASES.each do |test_entry|
    define_method("test_#{test_entry[:name]}") do
      trim = test_entry[:trim].nil? ? true : false
      raise_errors = test_entry[:raise_errors].nil? ? true : false
      result = CfeGotoh.frame_align(
        test_entry[:std],
        test_entry[:query],
        3,
        1,
        test_entry[:common_insert_locations],
        trim,
        raise_errors,
        true
      )
      assert_equal(test_entry[:expected_std], result[0])
      assert_equal(test_entry[:expected_query], result[1])
    end
  end

  ALIGNMENT_TEST_CASES = [
    {
      name: 'edges_are_trimmed',
      std: 'ACGTACGTACGT',
      query: 'CCCCCCACGTACGTACCTAAAAAA',
      expected_std: 'ACGTACGTACGT',
      expected_query: 'ACGTACGTACGT'
    },
    {
      name: 'edges_not_trimmed_when_not_specified',
      std: 'ACGTACGTACGT',
      query: 'AAAAAAACGTACGTACGTAAAAAA',
      trim: false,
      expected_std: '------ACGTACGTACGT------',
      expected_query: 'AAAAAAACGTACGTACGTAAAAAA'
    },
    {
      # std: ACGTGACGT-ACGT
      # qry: ACGT-ACGTGACGT
      name: 'indels_are_merged',
      std: 'ACGTGACGTACGT',
      query: 'ACGTACGTGACGT',
      expected_std: 'ACGTACGTACGT',
      expected_query: 'ACGTACGTACGT'
    },
    {
      # std: ACG--TA-CGTACGT
      # qry: ACGCCTATCGTACGT
      name: 'insertions_are_clustered',
      std: 'ACGTACGTACGT',
      query: 'ACGCCTATCGTACGT',
      expected_std: 'ACG---TACGTACGT',
      expected_query: 'ACGCCTATCGTACGT'
    },
    {
      # std: ACGCCTATCGTACGT
      # qry: ACG--TA-CGTACGT
      name: 'deletions_are_clustered',
      std: 'ACGCCTAGCGTACGT',
      query: 'ACGTACGTACGT',
      expected_std: 'ACGCCTAGCGTACGT',
      expected_query: 'ACG---TACGTACGT'
    },
    {
      # std: ACGT---ACGTACGT
      # qry: ACGTCCCACGTACGT
      name: 'insertions_are_frame_aligned',
      std: 'ACGTACGTACGT',
      query: 'ACGTCCCACGTACGT',
      expected_std: 'ACG---TACGTACGT',
      expected_query: 'ACGTCCCACGTACGT'
    },
    {
      # std: ACGTCCCACGTACGT
      # qry: ACGT---ACGTACGT
      name: 'deletions_are_frame_aligned',
      std: 'ACGTCCCACGTACGT',
      query: 'ACGTACGTACGT',
      expected_std: 'ACGTCCCACGTACGT',
      expected_query: 'ACG---TACGTACGT'
    },
    {
      # std: ACGTCCCACGTACGT
      # qry: ACGT---ACGTACGT
      name: 'insertions_moved_to_common_positions',
      std: 'ACGTCCCACGTACGT',
      query: 'ACGTACGTACGT',
      common_insert_locations: [3],
      expected_std: 'ACGTCCCACGTACGT',
      expected_query: 'ACGTACGTA---CGT',
    }
  ]

  ALIGNMENT_TEST_CASES.each do |test_entry|
    define_method("test_#{test_entry[:name]}") do
      trim = test_entry[:trim].nil? ? true : false
      raise_errors = test_entry[:raise_errors].nil? ? true : false
      result = CfeGotoh.frame_align(
        test_entry[:std],
        test_entry[:query],
        3,
        1,
        test_entry[:common_insert_locations],
        trim,
        raise_errors,
        false
      )
      assert_equal(test_entry[:expected_std], result[0])
      assert_equal(test_entry[:expected_query], result[1])
    end
  end

end
