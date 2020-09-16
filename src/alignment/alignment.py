from src.scoring import scoring
from src.objects import Spectrum, SequenceAlignment, HybridSequenceAlignment, Database, Alignments
from src.alignment import alignment_utils, hybrid_alignment
from src.sequence import gen_spectra

from src import utils
from src import database

import math
import re

####################### Time constants #######################

import time
TIME_LOG_FILE = './timelog.txt'

# to keep track of the time each step takes
FILTER_TIME = 0
FIRST_ALIGN_TIME = 0
AMBIGUOUS_REMOVAL_TIME = 0
PRECURSOR_MASS_TIME = 0
OBJECTIFY_TIME = 0

# keep track of frequency
FILTER_COUNT = 0
FIRST_ALIGN_COUNT = 0
AMBIGUOUS_REMOVAL_COUNT = 0
PRECURSOR_MASS_COUNT = 0
OBJECTIFY_COUNT = 0

####################### Public functions #######################

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
    
    # if EVERY position of seq1 is to the RIGHT of ALL positions of seq2, make a hybrid 
    if all([r < l for r in right_start for l in left_start]):
        return hybrid_alignment.hybrid_alignment(seq1, seq2)
    
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
        return hybrid_alignment.hybrid_alignment(seq1, seq2)
    
    # we have at least one good one. Return the shortest one
    nonhybrid_alignments.sort(key=lambda x: len(x))
    return (nonhybrid_alignments[0], None)

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

                    # append any alignments made from these 2 sequences
                    spec_alignments.append(
                        same_protein_alignment(b_seq, y_seq, prot_seq)
                    )
                
                # try just a dumb hybrid too to make sure
                spec_alignments.append((f'{b_seq}{y_seq}', f'{b_seq}-{y_seq}'))

            # otherwise try hybrid alignment
            else: 
                spec_alignments.append(hybrid_alignment.hybrid_alignment(b_seq, y_seq))
        
    # remove repeats
    return list(set([x for x in spec_alignments if x is not None]))

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
    global FIRST_ALIGN_TIME, AMBIGUOUS_REMOVAL_TIME, FILTER_TIME, PRECURSOR_MASS_TIME, OBJECTIFY_TIME
    global FIRST_ALIGN_COUNT, AMBIGUOUS_REMOVAL_COUNT, FILTER_COUNT, PRECURSOR_MASS_COUNT, OBJECTIFY_COUNT

    # run the first round of alignments
    st = time.time()
    a = align_b_y(b_hits, y_hits, db)

    FIRST_ALIGN_COUNT += len(a)
    FIRST_ALIGN_TIME += time.time() - st
    DEBUG and print(f'First alignment round took {time.time() - st} time resulting in {len(a)} alignments')

    # get the predicted length of the sequence and allow for a 25% gap to be filled in
    predicted_len = utils.predicted_len(max(spectrum.spectrum))
    allowed_gap = math.ceil(predicted_len * .25)

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
            alignment_utils.fill_in_precursor(spectrum, sequence, db, gap=allowed_gap, tolerance=precursor_tolerance) \
            if x is not None
        ]

        precursor_matches += p_ms

    PRECURSOR_MASS_COUNT += len(precursor_matches)
    PRECURSOR_MASS_TIME += time.time() - st
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
    updated_hybrids = [] if len(hyba) == 0 else hybrid_alignment.replace_ambiguous_hybrids(hyba, db, spectrum)

    AMBIGUOUS_REMOVAL_COUNT += len(updated_hybrids)
    AMBIGUOUS_REMOVAL_TIME += time.time() - st
    DEBUG and print(f'Getting rid of ambiguous time took {time.time() - st}')
 
    # Make alignments into the namedtuple types SpectrumAlignments
    # and HybridSequenceAlignments
    alignments = []
    tracker = {}
    st = time.time()
    for aligned_pair in nonhyba + updated_hybrids:

        # add to tracker, continue if what we see is already in the tracker
        if aligned_pair[0] in tracker:
            continue

        tracker[aligned_pair[0]] = True

        # get the precursor distance. If its too big, continue
        p_d = scoring.precursor_distance(spectrum.precursor_mass, gen_spectra.get_precursor(aligned_pair[0]))

        # get the final score of these sequences
        b_score = scoring.score_sequence(
            spectrum.spectrum, 
            sorted(gen_spectra.gen_spectrum(aligned_pair[0], ion='b')['spectrum']), 
            ppm_tolerance
        )
        y_score = scoring.score_sequence(
            spectrum.spectrum, 
            sorted(gen_spectra.gen_spectrum(aligned_pair[0], ion='y')['spectrum']), 
            ppm_tolerance
        )
        t_score = b_score + y_score

        # get the parent proteins of the sequence
        parents = alignment_utils.get_parents(aligned_pair[0] if aligned_pair[1] is None else aligned_pair[1], db)

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
    OBJECTIFY_COUNT += len(nonhyba + updated_hybrids)
    OBJECTIFY_TIME += time.time() - st
    DEBUG and print(f'Time to make into objects took {time.time() - st}')

    # print()
    # for a in sorted(alignments, key=lambda x: x.total_score, reverse=True)[:10]:
    #     print(f'{a.sequence} \t total: {a.total_score} \t b: {a.b_score} \t y: {a.y_score} \t hybrid: {"True" if len(a) == 8 else "False"}')

    if is_last:
        with open(TIME_LOG_FILE, 'w') as o:
            o.write(f'Total result filtering time: {FILTER_TIME}s \t seconds/op: {FILTER_TIME/FILTER_COUNT}s\n')
            o.write(f'B and Y alignment time: {FIRST_ALIGN_TIME}s \t seconds/op: {FIRST_ALIGN_TIME/FIRST_ALIGN_COUNT}s\n')
            o.write(f'Removing ambiguous hybrids time: {AMBIGUOUS_REMOVAL_TIME}s \t seconds/op: {AMBIGUOUS_REMOVAL_TIME/AMBIGUOUS_REMOVAL_COUNT}s\n')
            o.write(f'Matching precursor matches time: {PRECURSOR_MASS_TIME}s \t seconds/op: {PRECURSOR_MASS_TIME/PRECURSOR_MASS_COUNT}s\n')
            o.write(f'Turning matches into objects time: {OBJECTIFY_TIME} \t seconds/op: {OBJECTIFY_TIME/OBJECTIFY_COUNT}\n')

    return Alignments(
        spectrum, 
        sorted(
            alignments, 
            key=lambda x: (x.total_score, 1/x.precursor_distance, 1 if len(x) == 6 else 0), 
            reverse=True
        )[:n]
    )
