from gotoh import align_it, align_it_aa, align_it_aa_rb

# To test:
# * align_it_aa_rb
HXB2_NUCS = 'TGGAAGGGCTAATTCACTCCCAACGAAGACAAGATATCCTTGATCTGTGGATCTACCACACACAA'
DEFAULT_PENALTIES = [
    15,  # gap_open_penalty
    5,  # gap_extend_penalty
    1]  # use_terminal_gap_penalty


def test_align_it_perfect_match():
    hxb2 = HXB2_NUCS
    ref = seq = hxb2
    aligned_ref, aligned_seq, score = align_it(ref, seq, *DEFAULT_PENALTIES)
    assert aligned_ref == hxb2
    assert aligned_seq == hxb2
    assert score == 5 * len(hxb2)


def test_align_it_mismatch():
    ref = 'TGGAAGGGCTAATTCACTCCCAACGAAGACAAGATATCCTTGATCTGTGGATCTACCACACACAA'
    seq = 'TGGAAGGGATAATTCACTCCCAACGAAGACAAGATATCCTTGATCTGTGGATCTACCACACACAA'
    # mismatch:    ^
    expected_matches = len(ref) - 1
    expected_mismatches = 1
    aligned_ref, aligned_seq, score = align_it(ref, seq, *DEFAULT_PENALTIES)
    assert aligned_ref == ref
    assert aligned_seq == seq
    assert score == 5 * expected_matches - 4 * expected_mismatches


# noinspection DuplicatedCode
def test_align_it_insert():
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

    aligned_ref, aligned_seq, score = align_it(ref,
                                               seq,
                                               gap_open_penalty,
                                               gap_extend_penalty,
                                               use_terminal_gap_penalty)
    assert aligned_ref == expected_aref
    assert aligned_seq == seq
    assert score == (5 * expected_matches -
                     expected_gap_count * gap_open_penalty -
                     expected_gap_size * gap_extend_penalty)


# noinspection DuplicatedCode
def test_align_it_deletion():
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

    aligned_ref, aligned_seq, score = align_it(ref,
                                               seq,
                                               gap_open_penalty,
                                               gap_extend_penalty,
                                               use_terminal_gap_penalty)
    assert aligned_ref == ref
    assert aligned_seq == expected_aseq
    assert score == (5 * expected_matches -
                     expected_gap_count * gap_open_penalty -
                     expected_gap_size * gap_extend_penalty)


# noinspection DuplicatedCode
def test_align_it_aa_match():
    ref1 = 'R'
    seq1 = 'R'
    ref2 = 'D'
    seq2 = 'D'
    # align_it_aa uses different scores for each amino acid combination.
    expected_score1 = 7
    expected_score2 = 8

    aref1, aseq1, score1 = align_it_aa(ref1, seq1, *DEFAULT_PENALTIES)
    aref2, aseq2, score2 = align_it_aa(ref2, seq2, *DEFAULT_PENALTIES)

    assert (aref1, aseq1) == (ref1, seq1)
    assert (aref2, aseq2) == (ref2, seq2)
    assert score1 == expected_score1
    assert score2 == expected_score2


# noinspection DuplicatedCode
def test_align_it_aa_mismatch():
    ref0 = 'RRRRRRR'
    seq1 = 'RRRNRRR'
    seq2 = 'RRRDRRR'
    # align_it_aa uses different scores for each amino acid combination.
    expected_score1 = 6*7 - 5
    expected_score2 = 6*7 - 11

    aref1, aseq1, score1 = align_it_aa(ref0, seq1, *DEFAULT_PENALTIES)
    aref2, aseq2, score2 = align_it_aa(ref0, seq2, *DEFAULT_PENALTIES)

    assert (aref1, aseq1) == (ref0, seq1)
    assert (aref2, aseq2) == (ref0, seq2)
    assert score1 == expected_score1
    assert score2 == expected_score2


# noinspection DuplicatedCode
def test_align_it_aa_rb_mismatch():
    ref = 'RRRRRRR'
    seq = 'RRRDRRR'
    # Main difference is that align_it_aa_rb uses same score for all combos.

    aref, aseq = align_it_aa_rb(ref, seq, *DEFAULT_PENALTIES[:-1])

    assert aref == ref
    assert aseq == seq
