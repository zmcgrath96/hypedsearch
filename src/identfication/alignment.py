from src.utils import insort_by_index, all_perms_of_s, make_sparse_array
from src.scoring.scoring import score_subsequence, backbone_score, precursor_distance, xcorr, ion_backbone_score, intensity_ion_backbone_score, intensity_backbone_score, ion_intensity_percentage
from src.objects import Spectrum, SequenceAlignment, HybridSequenceAlignment, Database, Alignments
from src import database
from src.sequence.gen_spectra import get_precursor

from collections import defaultdict
import re 
from operator import itemgetter
from more_itertools import flatten
import math

import time

time_log_file = './timelog.txt'
filter_time = 0
first_align_time = 0
ambiguous_removal_time = 0
precursor_mass_time = 0
objectify_time = 0

filter_count = 0
first_align_count = 0
ambiguous_removal_count = 0
precursor_mass_count = 0
objectify_count = 0

hyb_alignment_pattern = re.compile(r'[-\(\)]')

########################################################################
#           TWO STRING ALIGNMENT FUNCTIONS
########################################################################

def align_overlaps(seq1: str, seq2: str) -> str:
    '''
    Attempt to align two string sequences. It will look at the right side of seq1 and left side of seq2
    to overlap the two strings. If no overlap is found, seq2 is appended to seq1

    Example 1: strings with overlap

        seq1: ABCD
        seq2: CDEF

        returns: ABCDEF

    Example 2: strings with no overlap

        seq1: ABCD
        seq2: EFGH

        return: ABCD-EFGH
    
    Inputs:
        seq1:    (str) the left side sequence 
        seq2:    (str) the right side sequence
    Outputs:
        str the attempted alignments
    '''
    alignment = None
    # if we have a perfect overlap, return it
    if seq1 == seq2:
        return seq1
    
    # if one is a full subsequence of another, return the larger one
    elif seq1 in seq2:
        return seq2
    elif seq2 in seq1:
        return seq1
    
    # try and find an alignment. seq2 should overlap as much of the right of seq1 as possible
    # get the starting points. 
    # Starting points means we found the first character in seq2 in seq1
    start_points = [i for i in range(len(seq1)) if seq2[0] == seq1[i]]
    for sp in start_points:

        # try and see if extending it makes it match
        # a correct overlap should mean we run out of characters in 
        # seq1 before we hit the end of seq2
        for i in range(sp, len(seq1)):

            # if the next value i in seq1 does not equal the next 
            # characeter i-sp in seq2
            if seq1[i] != seq2[i-sp]:
                i -= 1
                break

        # if i hits the the end of seq1, we have an overlap
        if i == len(seq1) - 1:
            s2_start = len(seq1) - sp
            right_seq = seq2[s2_start:] if s2_start < len(seq2) else ''
            alignment = seq1 + right_seq
            break
  
    # if no overlpa exists, just make append seq2 to seq1
    if alignment is None:
        alignment = seq1 + '-' + seq2

    return alignment

def same_protein_alignment(seq1: str, seq2: str, parent_sequence: str) -> (str, str):
    '''
    Attempt to create a non-hybrid alignment from two sequences from the same protein. If the 
    two sequences do not directly overlap but have <= 1% (or 1 if very short) of the total number 
    of amino acids, make the alignment. If not, create a hybrid alignment from that. If one 
    compeletely overlaps the other, use that as the alignment. 
    
    Example 1: Overlapped sequences
        seq1: CDE        starting position: 2
        seq2: ABCDEFG    starting position: 0
        
        Output: (ABCDEFG, None)
        
    Example 2: Partial overlaping sequences
        seq1: ABCDE      starting position: 0
        seq2: DEFGH      starting position: 3
        
        Outputs: (ABCDEFGH, None)
        
    Example 3: Non overlapping within the 1% (or 1)
        seq1: ABCDE      starting position: 0
        seq2: GHIJK      starting position: 6
        
        protein length: 52
        grace length = max(1, 1% of 52) = 1
        
        Outputs: (ABCDEFGHIJK, None)
        
    Example 4: Theoretical overlap but outside of range
        seq1: ABCDE      starting position: 0
        seq2: EFGHI      starting position: 30
        
        outside of allowed range, make hybrid
        
        Outputs: (ABCDEFGHI, ABCD(E)FGHI)
        
    
    Inputs:
        seq1:            (str) left sequence to align
        seq2:            (str) right sequence to align
        parent_sequence: (str) sequence of the shared parent protein
    Outputs:
        tuple:   first entry is the seqence, second entry is 
                 the second entry is the hybrid sequence, if not hybrid, then its None
    ''' 
    # check to see if they are equal or one covers the entirety of the other
    if seq1 == seq2:
        return (seq1, None)
    
    if seq1 in seq2:
        return (seq2, None)
    
    if seq2 in seq1: 
        return (seq1, None)
    
    # get the number of gap amino acids allowed
    gap_aa = max(1, len(parent_sequence) // 100)
    
    # get the positions of the left sequence from the protein
    left_start = [m.start() for m in re.finditer(seq1, parent_sequence)]
    
    # get the positions of the right sequence from the protein
    right_start = [m.start() for m in re.finditer(seq2, parent_sequence)]
    
    # if all of the right starts are the the left of the left starts, make a hybrid alignment
    if all([r < l for r in right_start for l in left_start]):
        return hybrid_alignment(seq1, seq2)
    
    # check to see if any of the points are within the gap_aa limit
    nonhybrid_alignments = []

    for l in left_start:
        for r in right_start:
            
            # if the right is to the left of left, continue
            if r < l:
                continue
                
            # if the right start - left start plus the length of the subsequence 
            # is less than the gap, just take the starting position of l and ending position 
            # of seq2 to make the full alignment
            if r - (l + len(seq1)) <= gap_aa:
                overlapped = parent_sequence[l: r + len(seq2)]
                nonhybrid_alignments.append(overlapped)
        
    # if no nonhybrids could be made, return a hybrid alignment
    if len(nonhybrid_alignments) == 0:
        return hybrid_alignment(seq1, seq2)
    
    # we have at least one good one. Return the shortest one
    nonhybrid_alignments.sort(key=lambda x: len(x))
    return (nonhybrid_alignments[0], None)
            

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
    attempted_overlap = align_overlaps(seq1, seq2)

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


def align_b_y(b_results: list, y_results: list, db: Database) -> list:
    '''
    Take 2 lists of sequences: one from the N terminus side (b_results) and 
    one from the C terminus side (y_results). Lookup the sequences in the database.
    If the two sequences are from the same protein, try and overlap the two strings.
    If there is an overlap, return it. In all other situations, a hybrid alignment is
    returned instead.

    Inputs:
        b_results:  (list of str) sequences found from b hits
        y_results:  (list of str) sequences fround from y hits
        db:         (Database) source of the sequences
    Outputs:
        (list) tuples of aligned sequences. First entry is the nonhybrid, second (if hybrid)
                has the hybrid characters -(). If not hybrid, it is None
    '''
    # try and create an alignment from each b and y sequence. Add all of the individual ones to start off with
    spec_alignments = [(seq, None) for seq in b_results + y_results]

    for b_seq in b_results:

        # get all the b proteins
        b_proteins = database.get_proteins_with_subsequence(db, b_seq)

        for y_seq in y_results:

            # ge the y proteins
            y_proteins = database.get_proteins_with_subsequence(db, y_seq)
            
            # the sequence is from the same protein, try and overlap it
            if any([x in y_proteins for x in b_proteins]):

                # get each protein they have in common
                shared_prots = [x for x in y_proteins if x in b_proteins]
                
                # try each of them 
                for sp in shared_prots:

                    # get the sequence from the entry for alignment
                    prot_seq = database.get_entry_by_name(db, sp).sequence

                    # try all same protein alignments
                    spec_alignments.append(
                        same_protein_alignment(b_seq, y_seq, prot_seq)
                    )

            # otherwise try hybrid alignment
            else: 
                spec_alignments.append(hybrid_alignment(b_seq, y_seq))
        
    # remove repeats
    return list(set([x for x in spec_alignments if x is not None]))

########################################################################
#          / TWO STRING ALIGNMENT FUNCTIONS
########################################################################

########################################################################
#           SINGLE STRING ALIGNMENT FUNCTIONS
########################################################################

def __get_surrounding_amino_acids(parent_sequence: str, sequence: str, count: int) -> list:
    '''
    Get the amino acids that surround a sequence. Return the (left, right) count number of amino acids

    Inputs:
        parent_sequence:    (str) the protein sequence to pull from
        sequence:           (str) the subsequence we are looking for
        count:              (int) the number of amino acids to get on each side
    Outputs:
        (list) pairs (tuples) of flanking amino acids for each occurance of sequence in parent_sequence
    '''
    # keep track of the pairs
    flanking_pairs = []

    # get all the positions of sequence in the parent sequence
    occurances = [m.start() for m in re.finditer(sequence, parent_sequence)]

    for o in occurances:
        flanking_pairs.append(
            (parent_sequence[max(0, o - count): o], 
            parent_sequence[o + len(sequence): min(len(parent_sequence), o + len(sequence) + count)])
        )

    return flanking_pairs

def __add_amino_acids(spectrum: Spectrum, sequence: str, db: Database, gap=3, tolerance=1) -> list:
    '''
    Try and add amino acids to get the closest precursor mass

    Inputs:
        spectrum:   (Spectrum) the observed precursor mass
        seqeunce:   (str) the attempted alignment
        db:         (Database) holds the protein sequences
    kwargs:
        gap:        (int) the number of additions allowed. Default=3
        tolerance:  (float) the mass (in Da) tolerance to accept in a precursor distance. Default=1
    Outputs:
        (list) the sequence(s) with the closest precursor mass
    '''
    filled_in  = []

    # get the parents of the sequence(s) and add or subtract amino acids
    parents = get_parents(sequence, db)

    # if its hybrid, we should fill it in in all possible ways
    if hyb_alignment_pattern.findall(sequence):

        # go through each set of parents
        for l_p in parents[0]:
            for r_p in parents[1]:

                # if the () are in the sequence, only add left and right of the sequence, don't
                # touch the overlap
                if '(' in sequence or ')' in sequence: 

                    # get the left and right proteins. 3rd entry (index 2) is the sequece
                    l_p_s = database.get_entry_by_name(db, l_p).sequence
                    r_p_s = database.get_entry_by_name(db, r_p).sequence

                    # separate the left and the right sequences
                    left_seq = sequence[:sequence.index(')')].replace('(', '')
                    right_seq = sequence[sequence.index('(')+1:].replace(')', '')

                    # get the sequence until ) and get the leftmost amino acids
                    left_aas = [x[0] for x in __get_surrounding_amino_acids(
                        l_p_s,
                        left_seq,
                        gap
                    )]

                    # get the sequence after ( and get the rightmost amino acids
                    right_aas = [x[0] for x in __get_surrounding_amino_acids(
                        r_p_s,
                        right_seq,
                        gap
                    )]

                    # go through all sets of left and right flanking amino acids
                    for l_aa in left_aas:
                        for r_aa in right_aas:

                            # go through all possible combinations of the flanking pairs
                            for i in range(gap + 1):
                                for j in range(gap - i + 1):

                                    # get the new sequence and its precursor mass. Add it to filled in list
                                    new_seq = l_aa[gap-i:] + sequence + r_aa[:j]
                                    new_prec = get_precursor(new_seq.replace('(', '').replace(')', ''))

                                    # get the precursor distance, and if it is close enough, keep it
                                    p_d = precursor_distance(spectrum.precursor_mass, new_prec)

                                    if p_d <= tolerance:
                                        filled_in.append(new_seq)

                # try all flanking of left right middle amino acids
                else:
                    left_seq = sequence.split('-')[0]
                    right_seq = sequence.split('-')[1]

                    left_flanking_pairs = __get_surrounding_amino_acids(
                        database.get_entry_by_name(db, l_p).sequence,
                        left_seq, 
                        gap
                    )

                    right_flanking_pairs = __get_surrounding_amino_acids(
                        database.get_entry_by_name(db, r_p).sequence,
                        right_seq, 
                        gap
                    )

                    # go through all the possible sets of left right middle sequences
                    for left_pair in left_flanking_pairs:
                        for right_pair in right_flanking_pairs:

                            # go through all lengths and contributions
                            for i in range(gap + 1):
                                for j in range(gap - i + 1):
                                    for k in range(max(gap - i - j + 1, 0)):
                                        for n in range(max(gap - i - j - k + 1, 0)):

                                            # leftmost addition
                                            ll = left_pair[0][gap-i:]
                                            
                                            # left center addition
                                            lc = left_pair[1][:j]

                                            # right center addition
                                            rc = right_pair[0][gap-k:]

                                            # right most addition
                                            rr = right_pair[1][:n]

                                            # get the new sequnce and precursor
                                            new_seq = ll + left_seq + lc + '-' + rc + right_seq + rr
                                            new_prec = get_precursor(new_seq.replace('-', ''))

                                            # get the precursor distance, and if it is close enough, keep it
                                            p_d = precursor_distance(spectrum.precursor_mass, new_prec)

                                            if p_d <= tolerance:
                                                                            filled_in.append(new_seq)
        

    # if its nonhybrid, try left and right side
    else:
        for p in parents[0]:

            # get the parent sequence
            p_seq = database.get_entry_by_name(db, p).sequence

            # get the flanking amino acid pairs
            for flanking_pair in __get_surrounding_amino_acids(p_seq, sequence, gap):

                # go through all possible combinations of the flanking pairs
                for i in range(gap + 1):
                    for j in range(gap - i + 1):

                        # get the new sequence and its precursor mass. Add it to filled in list
                        new_seq = flanking_pair[0][gap-i:] + sequence + flanking_pair[1][:j]

                        new_prec = get_precursor(new_seq)

                        # get the precursor distance, and if it is close enough, keep it
                        p_d = precursor_distance(spectrum.precursor_mass, new_prec)

                        if p_d <= tolerance:
                                    filled_in.append(new_seq)
    return filled_in

def __remove_amino_acids(spectrum: Spectrum, sequence: str, gap=3, tolerance=1) -> list:
    '''
    Remove up to gap number of amino acids to try and match precursor mass

    Inputs:
        spectrum:   (Spectrum) the aligned spectrum
        sequence:   (str) the attempted string alignment
    kwargs:
        gap:        (int) the total number of free amino acids to try. Default=3
        tolerance:  (float) the number (in Da) to allow as the tolerance for acceptable 
                            precursor masses. Default=1
    Outputs:
        (list) the sequence(s) with the closest precursor mass
    '''
    attempted = []

    # remove left and right first
    for i in range(gap + 1):
        for j in range(gap - i + 1):

             # if appended hybrid, removed the insides too
            if '-' in sequence:

                # split at the junction site
                left_seq = sequence.split('-')[0]
                right_seq = sequence.split('-')[1]

                for k in range(gap - i - j + 1):
                    for n in range(gap - i -j - k + 1):
                        
                        # create a new one by subtracting amino acids from each side
                        new_seq = left_seq[i:-j] if j >0 else left_seq [i:] \
                            + '-' + right_seq[k:-n] if n > 0 else right_seq[k:]

                        # find the new precursor mass
                        new_prec = get_precursor(new_seq.replace('-', ''))

                        # get the precursor distance, and if it falls within the tolerance, keep it
                        p_d = precursor_distance(spectrum.precursor_mass, new_prec)
                        
                        if p_d <= tolerance:
                                    attempted.append(new_seq)

            else:
                
                # take off the trailing (right) or leading (left) amino acids
                new_seq = sequence[i:-j] if j > 0 else sequence[i:]

                # if it s a hybrid with an ambiguous area, only subtract from the outside
                clean_seq = new_seq.replace('(', '').replace(')', '') \
                            if '(' in new_seq or ')' in new_seq else new_seq

                new_prec = get_precursor(clean_seq)

                # get the precursor distance, and if it is close enough, keep it
                p_d = precursor_distance(spectrum.precursor_mass, new_prec)

                if p_d <= tolerance:
                    attempted.append(new_seq)

    return attempted 



def fill_in_precursor(spectrum: Spectrum, sequence: str, db: Database, gap=3, tolerance=1) -> list:
    '''
    Try and fill in the gaps of an alignment. This is primarily focused on 
    filling in the gaps left by the difference in precursor mass. If we find that
    the difference is more than GAP amino acids, then return the broken one

    Inputs:
        spectrum:   (Spectrum) sequence to align 
        sequence:   (str) the aligned sequence
        db:         (Database) the database with string sequences
    kwargs:
        gap:        (int) the number of amino acids to accept. Default=3
        tolerance:  (float) the mass (in Da) to accept as error in matching precursors. Default=1
    Outputs:
        (list) sequence(s) that may have filled in the gap
    '''
    # the max mass for an amino acid (doubly charged)
    max_mass = 93.04

    # the min mass for an amino acid (doubly charged)
    min_mass = 28.5108

    # remove special characters for generating sequences
    clean_seq = sequence.replace('-', '').replace('(', '').replace(')', '')

    # get the theoretical precursor mass
    theory_precrusor = get_precursor(clean_seq)

    # get the precursor distance 
    p_d = precursor_distance(spectrum.precursor_mass, theory_precrusor)

    # if there are too many to add or subtract, return None
    if gap * max_mass < p_d:
        return [None]

    # if we can't add or subtract because our mass is too close, check to see if it
    # falls within our tolerance
    if p_d < min_mass:

        # if we are within the tolerance, return the sequence
        if p_d <= tolerance:
            return [sequence]

        # otherwise return none
        return [None]

    # determine HOW many amino acids we could be off
    num_off = min(math.ceil(p_d / max_mass) + 1, gap)

    # add amino acids
    if spectrum.precursor_mass > theory_precrusor:

        return __add_amino_acids(spectrum, sequence, db, num_off, tolerance)

    # subtract amino acids:
    else:

        return __remove_amino_acids(spectrum, sequence, num_off, tolerance)


def get_parents(seq: str, db: Database) -> (list, list):
    '''
    Get the parents of a sequence. If the sequence is a hybrid sequence, 
    then the second entry of the tuple holds a list of proteins for the right contributor.
    Otherwise the right entry is empty.

    Example 1: non hybrid peptide
        sequence: ABCDE
        Output: ([protein1, protein2], None)

    Example 2: hybridpeptide
        sequence: ABC(DE)FGH
        output: ([protein1], [protein2])

    Inputs:
        seq:    (str) sequence to find parents for
        db:     (Database) holds source information
    Outputs:
        (list, list) lists of parents
    '''
    get_sources = lambda s: database.get_proteins_with_subsequence(db, s)

    # If the sequence is hybrid, split it to find each parent
    if hyb_alignment_pattern.findall(seq):
        
        # straightforward left right split
        if '-' in seq:
            div = seq.split('-')
            left, right = div[0], div[1]
            return (get_sources(left), get_sources(right))

        # split by the () and make sure each side has the ambiguous area
        else:
            try:
                left = seq[:seq.index(')')].replace('(', '')
                right = seq[seq.index('(')+1:].replace(')', '')
                return (get_sources(left), get_sources(right))
            except:
                return ([], [])

    # otherwise just look up the one for the full sequence
    return (get_sources(seq), None)


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
    if len(database.get_proteins_with_subsequence(db, nonhyb)):
        return ((nonhyb, None))

    # Try replacing all L with I and vice versa
    possible = all_perms_of_s(nonhyb, 'LI')

    # if we had no other permutations, return the hybrid
    if len(possible) == 0:
        return hybrid

    # otherwise look through each permutation
    candidates = []

    for p in possible:

        # if this permutation exists as a nonhybrid, return it
        if len(database.get_proteins_with_subsequence(db, p)):
            return ((p, None))

        # otherwise keep it
        candidates.append(p)

    # take the highest scoring sequence
    best_candidates = sorted(
        candidates, 
        key=lambda x: sum(score_subsequence(observed.spectrum, p)), 
        reverse=True
    )

    # find the sequence that COULD exist
    for best_candidate in best_candidates:

        # find the hybrid junction
        for i in range(1, len(best_candidate)-1):
            left = best_candidate[:i]
            right = best_candidate[i:]

            # if we can find a left and a right with the sequence, return it
            if len(database.get_proteins_with_subsequence(db, left)) \
                and len(database.get_proteins_with_subsequence(db, right)):
                return ((best_candidate, f'{left}-{right}'))

    # worst case, return the hybrid
    return hybrid


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

########################################################################
#          / SINGLE STRING ALIGNMENT FUNCTIONS
########################################################################

def attempt_alignment(
    spectrum: Spectrum, 
    db: Database, 
    b_hits: list,
    y_hits: list, 
    n=3, 
    ppm_tolerance=20, 
    precursor_tolerance=1,
    scoring_alg='ibb', 
    DEBUG=False, 
    is_last=False
) -> list:
    '''
    Given a set of left and right (b and y ion) hits, try and overlap or extend one side to 
    explain the input spectrum

    Example:
        b_hits = [ABCD, LMNOP]
        y_hits = [ABCDEF]
        if the true value of the spectrum is ABCDEF, then we filter out LMNOP,
        overlap ABCD and ABCDEF and return ABCDEF

    Inputs:
        spectrum:               (Spectrum) spectrum to align
        db:                     (Database) Holds protein entries
        hits:                   (KmerMassesResults) hits from the hashing on a KmerMasses object
        base_kmer_len:          (int) minimum length kmer length used for filtering results
    kwargs:
        n:                      (int) number of results to return. Default=3
        ppm_tolerance:          (int) ppm tolerance to allow when scoring. Default=20
        precursor_tolerance:    (float) tolerance in Da to allow for a match. Default = 1
        scoring_alg:            (str) scoring algorithm to use. 'bb' for backbone, 'ion' for 
                              specific ion scoring, 'ibb' for ion backbone. Default='ibb'
    Outputs:
        attempted_alignments: (list) attempted alignemnts. Contains both or either of SequenceAlignment and HybridSequenceAlignment
    '''
    global first_align_time, ambiguous_removal_time, filter_time, precursor_mass_time, objectify_time
    global first_align_count, ambiguous_removal_count, filter_count, precursor_mass_count, objectify_count

    # align our b an y sequences
    st = time.time()
    a = align_b_y(b_hits, y_hits, db)

    first_align_count += len(a)
    first_align_time += time.time() - st
    DEBUG and print(f'First alignment round took {time.time() - st} time resulting in {len(a)} alignments')

    # Limit our search to things that match our precursor mass
    # try and fill in the gaps that are in any alignments
    st = time.time()
    precursor_matches = []

    for sequence_pairs in a:
        
        # take the sequence. If hybrid, take the hybrid, otherwise the non hybrid
        sequence = sequence_pairs[0] if sequence_pairs[1] is None else sequence_pairs[1]

        # add the closer precursors to the list
        p_ms = [
            x for x in \
            fill_in_precursor(spectrum, sequence, db, gap=4, tolerance=precursor_tolerance) \
            if x is not None
        ]

        precursor_matches += p_ms

    precursor_mass_count += len(precursor_matches)
    precursor_mass_time += time.time() - st
    DEBUG and print(f'Filling in precursor took {time.time() - st} for {len(a)} sequences')


    # seperate the hybrids from the non hybrids for later analysis
    nonhyba, hyba = [], []
    for p_m in precursor_matches:

        if '-' in p_m or '(' in p_m or ')' in p_m:
            hyba.append((p_m.replace('-', '').replace('(', '').replace(')', ''), p_m))
        else:
            nonhyba.append((p_m, None))
    
    # replace any hybrid alignments that are seen that can be explained by non 
    # hybrid sequences
    st = time.time()
    updated_hybrids = [] if len(hyba) == 0 else replace_ambiguous_hybrids(hyba, db, spectrum)

    ambiguous_removal_count += len(updated_hybrids)
    ambiguous_removal_time += time.time() - st
    DEBUG and print(f'Getting rid of ambiguous time took {time.time() - st}')
 
    # Make alignments into the namedtuple types SpectrumAlignments
    # and HybridSequenceAlignments
    alignments = []
    st = time.time()
    for aligned_pair in nonhyba + updated_hybrids:

        # get the precursor distance. If its too big, continue
        p_d = precursor_distance(spectrum.precursor_mass, get_precursor(aligned_pair[0]))

        # individual ion scores
        b_score = intensity_ion_backbone_score(spectrum, aligned_pair[0], 'b', ppm_tolerance)
        y_score = intensity_ion_backbone_score(spectrum, aligned_pair[0], 'y', ppm_tolerance)

        # backbone score
        bb_score = backbone_score(spectrum, aligned_pair[0], ppm_tolerance)

        # sum of the ion scores
        sion = sum(score_subsequence(spectrum.spectrum, aligned_pair[0], ppm_tolerance))

        # the ion intensity
        iip_score = ion_intensity_percentage(spectrum, aligned_pair[0], ppm_tolerance, '')


        t_score = bb_score // 8 + iip_score

        # backbone intensity
        bbi = intensity_backbone_score(spectrum, aligned_pair[0], ppm_tolerance)

        # get the parent proteins of the sequence
        parents = get_parents(aligned_pair[0] if aligned_pair[1] is None else aligned_pair[1], db)

        # check if the second entry is None
        if aligned_pair[1] is not None:
            alignments.append(
                HybridSequenceAlignment(
                    parents[0], 
                    parents[1], 
                    aligned_pair[0], 
                    aligned_pair[1], 
                    b_score, 
                    y_score, 
                    t_score, 
                    p_d
                )
            )

        # if its not a hybrid sequence, make a SequenceAlignment object
        else:
            alignments.append(
                SequenceAlignment(
                    parents[0], 
                    aligned_pair[0], 
                    b_score, 
                    y_score, 
                    t_score, 
                    p_d
                )
            )
    objectify_count += len(nonhyba + updated_hybrids)
    objectify_time += time.time() - st
    DEBUG and print(f'Time to make into objects took {time.time() - st}')

    tracker = {}
    set_alignments = []
    for a in alignments:
        if a.sequence in tracker: 
            continue
        tracker[a.sequence] = True
        set_alignments.append(a)
    # print(f'Number of alignments that passed the filter: {len(alignments)}')
    # print(sorted(alignments, key=lambda x: x.total_score, reverse=True)[:10])

    if is_last:
        with open(time_log_file, 'w') as o:
            o.write(f'Total result filtering time: {filter_time}s \t seconds/op: {filter_time/filter_count}s\n')
            o.write(f'B and Y alignment time: {first_align_time}s \t seconds/op: {first_align_time/first_align_count}s\n')
            o.write(f'Removing ambiguous hybrids time: {ambiguous_removal_time}s \t seconds/op: {ambiguous_removal_time/ambiguous_removal_count}s\n')
            o.write(f'Matching precursor matches time: {precursor_mass_time}s \t seconds/op: {precursor_mass_time/precursor_mass_count}s\n')
            o.write(f'Turning matches into objects time: {objectify_time} \t seconds/op: {objectify_time/objectify_count}\n')

    return Alignments(
        spectrum, 
        sorted(
            set_alignments, 
            key=lambda x: (x.total_score, x.b_score + x.y_score), 
            reverse=True
        )[:n]
    )
