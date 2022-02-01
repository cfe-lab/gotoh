import re
import gotoh
import webtypes


AlignItResult = webtypes.AlignItResult


class AlignIt :
    '''
    A class to call the Gotoh alignment library.

    Attributes
    ----------
    valid_nu : re.Pattern
        RegEx pattern to validate nucleotide sequences
    valid_aa : re.Pattern
        RegEx pattern to validate amino amino sequences

    Methods
    -------
    align_it(seqa, seqb, 
              gap_ini, gap_ext, 
              use_terminal_gap_penalty = False, 
              emulate_rb = False):
        Returns aligned sequences (with gaps) from the Gotoh algorithm.
        Expects nucleotide sequences.
    align_it_aa(seqa, seqb, 
                 gap_ini, gap_ext, 
                 use_terminal_gap_penalty = False, 
                 emulate_rb = False):
        Returns aligned sequences (with gaps) from the Gotoh algorithm.
        Expects amino acid sequences.
    '''
    def __init__(self):
        self.valid_nu = re.compile(r"^[acgturykmswbdhvnxACGTURYKMSWBDHVNX&]*$")
        self.valid_aa = re.compile(r"^[a-zA-Z&\*]*$")


    def align_it(self, seqa, seqb, 
                        gap_ini, gap_ext, 
                        use_terminal_gap_penalty = False, 
                        emulate_rb = False):
        '''
        Returns aligned sequences (with gaps) from the Gotoh algorithm.
        Expects nucleotide sequences, see align_it_aa() for amino acid.

                Parameters:
                        seqa (string): Nucleotide sequence (standard)
                        seqb (string): Another nucleotide sequence
                        gap_init (int): Gap initialization penalty
                        gap_extend (int): Gap extension penalty
                        use_terminal_gap_penalty (bool): penalize trailing gaps?
                        emulate_rb (bool): use original (Ruby) match/mismatch scores?

                Returns:
                        seqa (string): Aligned sequence a
                        seqb (string): Aligned sequence b
                        score (int): alignment score (gap penalties + match/mismatch)
                        exit_status (AlignItResult): ok, illegal_char, internal_error
        '''
        sa = ""
        sb = ""
        score = 0
        al_status = AlignItResult.internal_error

        try:
            if not bool(self.valid_nu.search(seqa) and self.valid_nu.search(seqb)):
                al_status = AlignItResult.illegal_char
            else:
                if emulate_rb:
                    [sa, sb] = gotoh.align_it_rb(seqa, seqb, gap_ini, gap_ext)
                    score = 0
                else:
                    [sa, sb, score] = gotoh.align_it(seqa, seqb, 
                                                     gap_ini, gap_ext, 
                                                     int(use_terminal_gap_penalty))
                al_status = AlignItResult.ok
        except:
            al_status = AlignItResult.internal_error

        return sa, sb, score, al_status


    def align_it_aa(self, seqa, seqb, 
                           gap_ini, gap_ext, 
                           use_terminal_gap_penalty = False,
                           emulate_rb = False):
        '''
        Returns aligned sequences (with gaps) from the Gotoh algorithm.
        Expects amino acid sequences, see align_it() for nucleotide.

                Parameters:
                        seqa (string): Amino acid sequence (standard)
                        seqb (string): Another amino acid sequence (seq)
                        gap_init (int): Gap initialization penalty
                        gap_extend (int): Gap extension penalty
                        use_terminal_gap_penalty (bool): penalize trailing gaps?
                        emulate_rb (bool): use original (Ruby) match/mismatch scores?

                Returns:
                        seqa (string): Aligned sequence a
                        seqb (string): Aligned sequence b
                        score (int): alignment score (gap penalties + match/mismatch)
                        exit_status (AlignItResult): ok, illegal_char, internal_error
        '''
        sa = ""
        sb = ""
        score = 0
        al_status = AlignItResult.internal_error

        try:
            if not bool(self.valid_nu.search(seqa) and self.valid_nu.search(seqb)):
                al_status = AlignItResult.illegal_char
            else:
                if emulate_rb:
                    [sa, sb] = gotoh.align_it_aa_rb(seqa, seqb, gap_ini, gap_ext)
                    score = 0
                else:
                    [sa, sb, score] = gotoh.align_it_aa(seqa, seqb, 
                                                     gap_ini, gap_ext, 
                                                     int(use_terminal_gap_penalty))
                al_status = AlignItResult.ok
        except:
            al_status = AlignItResult.internal_error

        return sa, sb, score, al_status