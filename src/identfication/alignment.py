from src.utils import insort_by_index, all_perms_of_s
from src.scoring.scoring import score_subsequence, backbone_score
from src.identfication.filtering import result_filtering
from src.types.objects import Spectrum, KmerMassesResults, SequenceAlignment, HybridSequenceAlignment
from src.types.database import Database

from collections import defaultdict
import re 

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
        alignment = seq1
    
    # if one is a full subsequence of another, return the larger one
    elif seq1 in seq2:
        alignment = seq2
    elif seq2 in seq1:
        alignment = seq1
    
    else:
        # try and find an alignment. seq2 should overlap as much of the right of seq1 as possible
        start_points = [i for i in range(len(seq1)) if seq2[0] == seq1[i]]
        for sp in start_points:
            # try and see if extending it makes it match
            for i in range(sp, len(seq1)):
                if seq1[i] != seq2[i-sp]:
                    i -= 1
                    break
            if i == len(seq1) - 1:
                s2_start = len(seq1) - sp
                right_seq = seq2[s2_start:] if s2_start < len(seq2) else ''
                alignment = seq1 + right_seq
                break
  
    if alignment is None:
        alignment = seq1 + '-' + seq2
    return alignment
            

def hybrid_alignment(seq1: str, seq2: str) -> (str, str):
    '''
    Create a hybrid alignment from 2 sequences. If an overlap between these two sequences
    is found, () are placed around the ambiguous section. If there is no overlap, then 
    seq2 is appended to seq1 with a - at the junction

    Example 1: overlap of two strings 

        seq1: ABCDE
        seq2: DEFGH

        attempted alignment: ABC(DE)FGH

    Example 2: no overlap between the two strings

        seq1: ABCD
        seq2: EFGH

        attempted alignment: ABCD-EFGH
        
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
    if attempted_overlap is not None and 0 < len(attempted_overlap) < len(seq1) + len(seq2):
        # there is an overlap and some ambiguity
        # get the starting point of seq2
        rightstart = attempted_overlap.index(seq2)
        leftend = len(seq1) - 1
        # range between leftend and rigth start is ambiguous
        middle_sec = attempted_overlap[rightstart:leftend + 1]
        
        alignment = attempted_overlap
        hybalignment = attempted_overlap[:rightstart] + '(' + middle_sec + ')' + attempted_overlap[leftend+1:]
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
        (list) matched b and y hits that form good alignments
    '''
    # try and create an alignment
    spec_alignments = []
    for bs in b_results:
        bproteins = [id_ for id_, _ in db.tree.find_all(bs)]
        for ys in y_results:
            yproteins = [id_ for id_, _ in db.tree.find_all(ys)]
            
            # the sequence is from the same protein, try and overlap it
            if any([x in yproteins for x in bproteins]):
                spec_alignments.append(align_overlaps(bs, ys))
            # otherwise try hybrid alignment
            else: 
                spec_alignments.append(hybrid_alignment(bs, ys)[1])
        
    # remove repeats
    return list(set([x for x in spec_alignments if x is not None]))

########################################################################
#          / TWO STRING ALIGNMENT FUNCTIONS
########################################################################

########################################################################
#           SINGLE STRING ALIGNMENT FUNCTIONS
########################################################################

def get_parents(seq: str, db: Database) -> (list, list):
    '''
    Get the parents of a sequence. If the sequence is a hybrid sequence, 
    then the second entry of the tuple holds a list of proteins for the right contributor.
    Otherwise the right entry is empty.

    Inputs:
        seq:    (str) sequence to find parents for
        db:     (Database) holds source information
    Outputs:
        (list, list) lists of parents
    '''
    get_sources = lambda s: [x[0] for x in db.tree.find_all(s)]

    if hyb_alignment_pattern.findall(seq):
        if '-' in seq:
            div = seq.split('-')
            left, right = div[0], div[1]
            return (get_sources(left), get_sources(right))
        else:
            left = seq[:seq.index(')')].replace('(', '')
            right = seq[seq.index('(')+1:].replace(')', '')
            return (get_sources(left), get_sources(right))
    else:
        return (get_sources(seq), None)


def replace_ambiguous_hybrids(hybrid_alignments: list, db: Database) -> list:
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
        hybrid_alignments:  (list of str) attempted hybrid alignments
        db:                 (Database) the source of the sequences
    Ouputs:
        (list) alignments. If no replacements are found, the output is the input 
    '''
    ret = []
    for hybalignment in hybrid_alignments:
        added = False
        nonhyb = re.sub(hyb_alignment_pattern, '', hybalignment)

        # also try to replace L and I with eachother
        possible = all_perms_of_s(nonhyb, 'LI')
        if len(possible) == 0 and db.tree.find(nonhyb):
            ret.append(nonhyb)
            added = True
        else:
            for p in possible:
                if db.tree.find(p):
                    ret.append(p)
                    added = True
                    break
            if not added:
                ret.append(hybalignment)

    return ret 

########################################################################
#          / SINGLE STRING ALIGNMENT FUNCTIONS
########################################################################

def attempt_alignment(spectrum: Spectrum, db: Database, hits: KmerMassesResults, base_kmer_len: int, n=3, ppm_tolerance=20, scoring_alg='ibb') -> list:
    '''
    Given a set of hits in the form of KmerMassesResult (lists of sequences), reduce
    the number of hits and filter out poor scoring sequences. Combine N and C terminus
    amino acid sequences (via string overlapping) to explain the spectrum.

    Example:
        hits = [ABCD, LMNOP, ABCDEF]
        if the true value of the spectrum is ABCDEF, then we filter out LMNOP,
        overlap ABCD and ABCDEF and return ABCDEF

    Inputs:
        spectrum:       (Spectrum) spectrum to align
        db:             (Database) Holds the tree and protein entries
        hits:           (KmerMassesResults) hits from the hashing on a KmerMasses object
        base_kmer_len:  (int) minimum length kmer length used for filtering results
    kwargs:
        n:              (int) number of results to return. Default=3
        ppm_tolerance:  (int) ppm tolerance to allow when scoring. Default=20
        scoring_alg:    (str) scoring algorithm to use. 'bb' for backbone, 'ion' for 
                              specific ion scoring, 'ibb' for ion backbone. Default='ibb'
    Outputs:
        attempted_alignments: (list) attempted alignemnts. Contains both or either of SequenceAlignment and HybridSequenceAlignment
    '''

    # narrow down the potential alignments
    b_results, y_results = result_filtering(spectrum, hits, base_kmer_len, scoring_alg=scoring_alg)

    # align our b an y sequences
    a = align_b_y(b_results, y_results, db)

    # sort them by their score against the spectrum and save them
    a = sorted(a, key=lambda x: backbone_score(spectrum, x.replace('-', '').replace('(', '').replace(')', ''), ppm_tolerance), reverse=True)
    
    # seperate the hybrids from the non hybrids for later analysis
    nonhyba, hyba = [], []
    for aligned in a:
        if hyb_alignment_pattern.findall(aligned):
            hyba.append(aligned)
        else:
            nonhyba.append(aligned)
    
    # replace any hybrid alignments that are seen that can be explained by non 
    # hybrid sequences
    updated_hybrids = [] if len(hyba) == 0 else replace_ambiguous_hybrids(hyba, db)

    # Make alignments into the namedtuple types SpectrumAlignments
    # and HybridSequenceAlignments
    alignments = []
    for sequence in nonhyba + updated_hybrids:
        parents = get_parents(sequence, db)
        if hyb_alignment_pattern.findall(sequence):
            aa_seq = re.sub(hyb_alignment_pattern, '', sequence)
            b_score, y_score = score_subsequence(spectrum.spectrum, aa_seq)
            bb_score = backbone_score(spectrum, aa_seq, ppm_tolerance)
            alignments.append(HybridSequenceAlignment(parents[0], parents[1], aa_seq, sequence, b_score, y_score, bb_score))
        else:
            b_score, y_score = score_subsequence(spectrum.spectrum, sequence)
            bb_score = backbone_score(spectrum, sequence, ppm_tolerance)
            alignments.append(SequenceAlignment(parents[0], sequence, b_score, y_score, bb_score))

    return sorted(alignments, key=lambda x: x.total_score, reverse=True)[:n]
