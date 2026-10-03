def global_alignment(seq1, seq2, scoring_function):
    """Global sequence alignment using the Needleman–Wunsch algorithm.

    Indels should be denoted with the "-" character.

    Parameters
    ----------
    seq1: str
        First sequence to be aligned.
    seq2: str
        Second sequence to be aligned.
    scoring_function: Callable

    Returns
    -------
    str
        First aligned sequence.
    str
        Second aligned sequence.
    float
        Final score of the alignment.

    Examples
    --------
    >>> global_alignment("abracadabra", "dabarakadara", lambda x, y: [-1, 1][x == y])
    ('-ab-racadabra', 'dabarakada-ra', 5.0)

    Other alignments are not possible.

    """
    seq1_length = len(seq1)
    seq2_length = len(seq2)

    # Initialise a matrix
    dp_matrix = []

    for i in range(seq1_length + 1):
        row_vals = []
        for j in range(seq2_length + 1):
            if i == 0:
                row_vals.append(j*-1)
            else:
                if j == 0:
                    row_vals.append(i*-1)
                else:
                    row_vals.append(0)
                
        dp_matrix.append(row_vals)

    for i in dp_matrix:
        print(i)

    gap_penalty = -1

    # Go over each element in the row and score it like 
    for i in range(1, seq1_length + 1):
        char_a = seq1[i - 1]
        for j in range(1, seq2_length + 1):
            char_b = seq2[j - 1]

            # Compare both of these using the scoring function - diagonal movement
            match_mistmatch_score = dp_matrix[i - 1][j - 1] + scoring_function(char_a, char_b)

            # gap from down
            gap_down = dp_matrix[i - 1][j] + gap_penalty

            # gap from left
            gap_left = dp_matrix[i][j - 1] + gap_penalty

            # Pick out the max value
            dp_matrix[i][j] = max(match_mistmatch_score, gap_down, gap_left)

    for i in dp_matrix:
        print(i)

    # start at the bottom corner and get a path
    seq_1_align = []
    seq_2_align = []

    i = seq1_length
    j = seq2_length

    final_score = 0

    while i > 0 or j > 0:
        score_max = 0

        # determine if gap from top is the best
        gap_top = dp_matrix[i - 1][j] + gap_penalty

        # determine if gap from left is the best
        gap_left = dp_matrix[i][j - 1] + gap_penalty

        # determine if the match/mismatch is the best
        match_mismatch = dp_matrix[i - 1][j - 1] + scoring_function(seq1[i - 1], seq2[j - 1])

        # Now get the score_max
        score_max = max(gap_left, gap_top, match_mismatch)

        # Based on the score, determine where to go
        
        if score_max == gap_top:
            # Go above (seq1 will change value but seq2 will not)
            i = i - 1
            seq_1_align.append(seq1[i])
            seq_2_align.append('-')
            final_score += gap_penalty

        elif score_max == gap_left:
            # Go left (seq2 will change vlaue but seq1 will not)
            j = j - 1
            seq_2_align.append(seq2[j])
            seq_1_align.append('-')
            final_score += gap_penalty
        else:
            # Diagnoal (both of them get a val)
            i = i - 1
            j = j - 1
            seq_1_align.append(seq1[i])
            seq_2_align.append(seq2[j])
            final_score += scoring_function(seq1[i], seq2[j])

    seq_1_align.reverse()
    seq_2_align.reverse()

    seq_1_align_res = "".join(seq_1_align)
    seq_2_align_res = "".join(seq_2_align)
    
    print(seq_1_align_res)
    print(seq_2_align_res)
    print(final_score)

global_alignment("abracadabra", "dabarakadara", lambda x, y: [-1, 1][x == y])

def local_alignment(seq1, seq2, scoring_function):
    """Local sequence alignment using the Smith-Waterman algorithm.

    Indels should be denoted with the "-" character.

    Parameters
    ----------
    seq1: str
        First sequence to be aligned.
    seq2: str
        Second sequence to be aligned.
    scoring_function: Callable

    Returns
    -------
    str
        First aligned sequence.
    str
        Second aligned sequence.
    float
        Final score of the alignment.

    Examples
    --------
    >>> local_alignment("pending itch", "unending glitch", lambda x, y: [-1, 1][x == y])
    ('ending --itch', 'ending glitch', 9.0)

    Other alignments are not possible.

    """
    raise NotImplementedError()


## This is an example scoring function, you should implement a version which uses a scoring matrix 
def scoring_function_simple(aa_i,aa_j):
    score = [-1, 1][aa_i == aa_j]
    return (score)
