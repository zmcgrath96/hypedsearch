from src.objects import Spectrum
from src.database import Database
from src.utils import all_perms_of_s
from src.scoring import scoring
from src.alignment import alignment_utils

from src import database

#################### Private functions ####################

def __replace_ambiguous_hybrid(hybrid: tuple, db: Database, observed: Spectrum) -> (str, str):
    '''
    Attempt to replace a hybrid with a sequence from the database. 

    Inputs:
        hybrid:     (tuple) (nonhybrid sequence, hybrid sequence)
        db:         (Database) the source of sequences
        observed:   (Spectrum) observed spectrum
    Outputs:
        (str, str) (the updated non hybrid sequence, None if not hybrid else hybrid sequence)
    '''
    # get the sequence without the hybrid characters -()
    nonhyb = hybrid[0]

    # see if the sequence exists as a non hybrid. If so, return that
    if len(db.get_proteins_with_subsequence(nonhyb)):
        return ((nonhyb, None))

    # Try replacing all L with I and vice versa
    possible = all_perms_of_s(nonhyb, 'LI')

    # if we had no other permutations, return the hybrid
    if len(possible) == 0:
        return hybrid

    # try and find a sequnce that could be explained by an LI switch
    for p in possible:

        # if this permutation exists as a nonhybrid, return it
        if len(db.get_proteins_with_subsequence(p)):
            return ((p, None))

    # if we didn't find a sequence that could be found in the database, 
    # just return the hybrid input 
    return hybrid

#################### Public functions ####################

def replace_ambiguous_hybrids(hybrid_alignments: list, db: Database, observed: Spectrum) -> list:
    '''
    Remove any ambiguous hybrid alignments that can be explained by non hybrid sequences.
    The returned list has the sequences or their replacements in the same order that 
    they were in on entry. 

    Amino acids L and I are swapped and tried in the search due to the ambiguity 
    in their mass

    Example:
        hybrid_alignments: [ABC(DE)FG, LMN-OPQ]

        ABCDEFG is found in the database as a non hybrid sequence, LMNOPQ is not

        returned list is then [ABCDEFG, LMN-OPQ]

    Inputs:
        hybrid_alignments:  (list) tuples of attempted hybrid alignments with (nonhyb, hyb) form
        db:                 (Database) the source of the sequences
        observed:           (Spectrum) the observed spectrum
    Ouputs:
        (list) alignments. If no replacements are found, the output is the input 
    '''
    return [
        __replace_ambiguous_hybrid(hybrid_alignment, db, observed)\
         for hybrid_alignment in hybrid_alignments
    ]

def hybrid_alignment(seq1: str, seq2: str) -> (str, str):
    '''
    Create a hybrid alignment from 2 sequences. If an overlap between these two sequences
    is found, () are placed around the ambiguous section. If there is no overlap, then 
    seq2 is appended to seq1 with a - at the junction

    Example 1: overlap of two strings 

        seq1: ABCDE
        seq2: DEFGH

        attempted alignment: ABC(DE)FGH
        Output: (ABCDEFGH, ABC(DE)FGH)

    Example 2: no overlap between the two strings

        seq1: ABCD
        seq2: EFGH

        attempted alignment: ABCD-EFGH
        Output: (ABCDEFGH, ABCD-EFGH)
        
    Inputs:
        seq1:    (str) the left sequence
        seq2:    (str) the right sequence
    Outputs:
        (str, str)  a tuple of strings. The first string is the sequence without any hybrid  
                    identifying strings. The second is the sequence with the hybrid
                    identifying strings -()
    '''
    alignment = ''
    hybalignment = ''
    attempted_overlap = alignment_utils.align_overlaps(seq1, seq2)

    # If an alignment was made (the shorter length means combined sequences is shorter than
    # just appending them), identify the ambiguous area
    if attempted_overlap is not None and 0 < len(attempted_overlap) < len(seq1) + len(seq2):
        
        # there is an overlap and some ambiguity
        # get the starting point of seq2
        rightstart = attempted_overlap.index(seq2)
        leftend = len(seq1) - 1

        # range between leftend and rigth start is ambiguous
        middle_sec = attempted_overlap[rightstart:leftend + 1]
        alignment = attempted_overlap
        hybalignment = attempted_overlap[:rightstart] \
                        + '(' + middle_sec + ')' \
                        + attempted_overlap[leftend+1:]
    
    # therwise just append right seq to the left with a - to seperate them 
    else:
        alignment = seq1 + seq2
        hybalignment = seq1 + '-' + seq2
        
    return (alignment, hybalignment)