from src.types.objects import KmerMassesResults, Spectrum
from src.scoring.scoring import score_subsequence, backbone_score, ion_backbone_score
from src.utils import insort_by_index

from statistics import mean
from typing import Iterable
from collections import defaultdict

def slope_filtering(a: Iterable, min_window_size=5, mean_filter=1, key=None) -> list:
    '''
    Filter out values by slope. Use a window size of at least min_window_size and at max 1% the size of 
    iterable a. Result is the subset of results whose window mean > mean_filter*mean_slope. a should be sorted
    with the steepest slope at the earliest indices. Key is used if provided to index, otherwise it is 
    assumed that the list elements are the scores

    Example:
        a: [10, 9, 8, 7, 6.3, 5.8, 5.4, 5.3, 5.25, 5.23]
        min window size: 1
        mean_filter: 1

        slopes of a = [1, 1, 1, .7, .5, .4, .1, .05, .02]
        mean slope of a = .53

        output is then [10, 9, 8, 7, 6.3]
    
    Inputs:
        a:     (Iterable) the list or list like object to filter
    kwargs:
        min_window_size:    (int) the minimum window size to use for sliding means. Default=5
        mean_filter:        (int) a scaling factor to determine what values pass. This value
                                  is mutliplied by the mean slope and values > than this are accepted. Default=1
    Outputs:
        list with the filtered results
    '''
    # we need data points. If not enough given, just return the list
    if len(a) < 2 * min_window_size:
        return a

    # make a minimum window size (for large lists) of 1%
    adjusted_window_size = len(a) // 100

    # try and make the window size about 1% of the size for large lists
    window_size = min_window_size if adjusted_window_size < min_window_size else adjusted_window_size
    
    # lambda to grab the score value
    get_val = lambda x: x if key is None else (x[key] if type(x) == dict or type(x) == list or type(x) == tuple else getattr(x, key))

    # calculate the slope for a window size
    slopes = [(get_val(a[i+window_size]) - get_val(a[i]))/window_size for i in range(len(a) - window_size)]
    avgslope = mean(slopes)
    filtered = []

    # keep anything with a slope > filter
    for i, slope in enumerate(slopes):

        # if the slope passes the filter, add those values
        if abs(slope) > abs(avgslope):
            filtered.append(a[i])
        
        # if the slope doesnt pass the filter, stop adding values
        else:
            break

    return filtered

def hash_and_score_base_kmers(
    spectrum: Spectrum, 
    hits: list, 
    base_kmer_length: int, 
    scoring_alg: callable, 
    ion: str, 
    ppm_tolerance: int
) -> (defaultdict, list):
    '''
    Take a set of initial hits, reduce the hits to lists based on their base kmers, and score the base kmers.
    Any base k-mer that does not get score > 0 is not added to the list
    Example:
        hits: [MAL, MALWAR, PPST, PPR]
        ion: b
        base_kmer_length: 3

        base kmer reduction: {MAL: [MAL, MALWAR], PPS: [PPST], PPR: [PPR]}
        scored base kmers: [(MAL, 1.1), (PPS, 0.8), (PPR, 0.7)]

        return value: ({MAL: [MAL, MALWAR], PPS: [PPST], PPR: [PPR]}, [(MAL, 1.1), (PPS, 0.8), (PPR, 0.7)])

    Ion type determines if base kmers are left to right or right to left. Above example uses b ion, so kmers are
    built left to right.

    Inputs:
        spectrum:           (Spectrum) used for scoring the base sequences
        hits:               (list) string sequences from initial hits
        base_kmer_length:   (int) the length of the base kmer to use for indexing
        scoring_alg         (callable) function used for scoring. Inputs should be, in order:
                                    spectrum: Spectrum, reference_sequence: str, ion: str, ppm_tolerance: int
        ion:                (str) the ion list to use. Options are [b, y]
        ppm_tolerance:      (int) the ppm_tolerance to accept while scoring
    Outputs:
        (defaultdict, list) basekmer hashed sequences, sequence score tuples list
    '''
    # keep track of bad sequences
    blacklist = {}

    # the binned sequences
    binned = defaultdict(list)

    # the list of sequence, score tuples
    scores = []

    # get the base mer based on ion
    get_base_mer = lambda seq: seq[:base_kmer_length] if ion == 'b' else seq[len(seq)-base_kmer_length:]

    for masssequence_hit in hits:
        base_mer = get_base_mer(masssequence_hit)

        # check to see if we've seen it
        # skip if its blacklisted
        if base_mer in blacklist:
            continue

        # if we've not seen it, try and add it
        if base_mer not in binned:
            base_mer_score = scoring_alg(spectrum, base_mer, ion, ppm_tolerance)

            # blacklist for bad scores
            if not base_mer_score > 0:
                blacklist[base_mer] = None
                continue

            # insert in order
            scores = insort_by_index((base_mer, base_mer_score), scores, 1)
        binned[base_mer].append(masssequence_hit)

    return (binned, scores)

def result_filtering(
    spectrum: Spectrum, 
    hits: KmerMassesResults, 
    base_kmer_length: int, 
    ppm_tolerance=20, 
    scoring_alg='ibb'
) -> (list, list):
    '''
    Filter out the sequences that do not score well. Use a mean filter to do so.
    Starting with KmerMassesResults namedtuple, each entry a list of MassSequence namedtuples,
    filter out by the base kmer length. This bins sequences 'MALWAR' into the same one as 'MALWARMQAAR'
    The point of doing so is to reduce the interested sequences. Once we filter out by the base kmers, 
    we then score all the remaining sequences and return the top scoring ones.

    Inputs:
        spectrum:           (Spectrum) the mass spectrum being aligned
        hits:               (KmerMassesResults) raw results from the search algorithm on the KmerMasses 
        base_kmer_length:   (int) the smallest length k-mer to consider
    kwargs:
        ppm_tolerance:      (int) the ppm tolerance to accept when scoring. Default=20
        scoring_alg:        (str) scoring algoirithm to use. 'bb' for backbone, 'ion' for separate ion score, 
                                  'ibb' for ion backbone. Default='ibb'
    Outputs:
        (list, list)    list of strings of the results from the left (N terminus) and right (C terminus) sides 
                        respectively
    '''
    # keep track of the sequences that start with the base kmers 
    base_mer_hashed_b = defaultdict(list)
    base_mer_hashed_y = defaultdict(list)
    
    # scoring algorithm to use
    def score_alg(spectrum, refseq, ion, ppm_tolerance):

        # ion backbone score
        if 'ibb' == scoring_alg:
            return ion_backbone_score(spectrum, refseq, ion, ppm_tolerance)

        # ion counting score
        if 'ion' == scoring_alg:
            retindex = 0 if ion == 'b' else 1
            return score_subsequence(spectrum.spectrum, refseq, ppm_tolerance=ppm_tolerance)[retindex]
        
        # default to backbone score
        return backbone_score(spectrum, refseq, ppm_tolerance)

    # sorted lists of (basekmers, score) by score
    b_scores = []
    y_scores = []

    # hash by the base kmer
    for hittype, hitlist in hits._asdict().items():
        if 'b' in hittype:

            new_binned_b, new_b_scores = hash_and_score_base_kmers(
                spectrum, 
                hitlist, 
                base_kmer_length, 
                score_alg,
                'b', 
                ppm_tolerance
            )
            
            # add these to the running values
            for k, v in new_binned_b.items():
                base_mer_hashed_b[k] += v

            b_scores += new_b_scores

        else:
            new_binned_y, new_y_scores = hash_and_score_base_kmers(
                spectrum, 
                hitlist, 
                base_kmer_length, 
                score_alg,
                'y', 
                ppm_tolerance
            )
            
            print('new binned y')
            print(new_binned_y)
            # add these to the running values
            for k, v in new_binned_y.items():
                base_mer_hashed_y[k] += v

            y_scores += new_y_scores

    toscoreb = []
    toscorey = []

    # reverse the sorting of score to go from high to low
    b_scores = b_scores[::-1]
    y_scores = y_scores[::-1]
    
    # mean filtered base kmers by scores
    b_filtered = slope_filtering(b_scores, key=1)
    y_filtered = slope_filtering(y_scores, key=1)
    
    if len(b_filtered) < 5: # default to be able to report some result
        b_filtered = b_scores[::-1][:5]
    if len(y_filtered) < 5:
        y_filtered = y_scores[::-1][:5]
    
    # revese to go from highest to lowest
    for basemerb in b_filtered:
        toscoreb += [basemerb[0]] + base_mer_hashed_b[basemerb[0]]
        
    for basemery in y_filtered:
        toscorey += [basemery[0]] + base_mer_hashed_y[basemery[0]]
    
    # sorted score of all leftover sequences from highest to lowest
    best_b_results = sorted(toscoreb, key=lambda mer: score_alg(spectrum, mer, 'b', ppm_tolerance), reverse=True)
    best_y_results = sorted(toscorey, key=lambda mer: score_alg(spectrum, mer, 'y', ppm_tolerance), reverse=True)

    return (best_b_results, best_y_results)