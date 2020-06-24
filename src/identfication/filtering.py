from src.types.objects import KmerMassesResults, Spectrum
from src.scoring.scoring import score_subsequence, backbone_score, ion_backbone_score, xcorr, overlap_score
from src.utils import insort_by_index, make_sparse_array
from src.spectra.gen_spectra import gen_spectrum

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

def mean_filtering(a: Iterable, mean_filter=1, key=None) -> list:
    '''
    Filter out values < the mean value * the mean filter

    Example:
        a: [10, 9, 8, 7, 6, 5]
        average: 7.5
        mean_filter: 1

        output: [10, 9, 8]

    Inputs:
        a:      (Iterable) values to filter
    kwargs:
        mean_filter:    (float) the number of means to accept by our filter. Default=1
        key:            (any) the way to index a value in the iterable if each entry is an object. Default=None
    Outputs:
        (list) values that pass the filter
    '''
    if len(a) == 0:
        return a

    # get the average value
    avg = mean([x[key] for x in a]) if key is not None else mean(a)

    # return anything above the filter
    return [x for x in a if x[key] >= mean_filter * avg] \
        if key is not None else [x for x in a if x >= mean_filter * avg]

def hash_base_kmers(
    spectrum: Spectrum, 
    hits: list, 
    base_kmer_length: int, 
    ion: str, 
) -> defaultdict:
    '''
    Take a set of initial hits, reduce the hits to lists based on their base kmers.
    Example:
        hits: [MAL, MALWAR, PPST, PPR]
        ion: b
        base_kmer_length: 3

        base kmer reduction: {MAL: [MAL, MALWAR], PPS: [PPST], PPR: [PPR]}

        return value: ({MAL: [MAL, MALWAR], PPS: [PPST], PPR: [PPR]}

    Ion type determines if base kmers are left to right or right to left. Above example uses b ion, so kmers are
    built left to right.

    Inputs:
        spectrum:           (Spectrum) used for scoring the base sequences
        hits:               (list) string sequences from initial hits
        base_kmer_length:   (int) the length of the base kmer to use for indexing
        ion:                (str) the ion list to use. Options are [b, y]
    Outputs:
        (defaultdict) basekmer hashed sequences,
    '''
    # keep track of bad sequences
    blacklist = {}

    # the binned sequences
    binned = defaultdict(list)

    # get the base mer based on ion
    get_base_mer = lambda seq: seq[:base_kmer_length] if ion == 'b' else seq[len(seq)-base_kmer_length:]

    for masssequence_hit in hits:
        base_mer = get_base_mer(masssequence_hit)

        # check to see if we've seen it
        # skip if its blacklisted
        if base_mer in blacklist:
            continue

        binned[base_mer].append(masssequence_hit)

    return binned

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

    # make a sparse spectrum array in case xcorr score is used
    sparse_spectrum = make_sparse_array(spectrum.spectrum, .02)
    
    # scoring algorithm to use
    def score_alg(spectrum, refseq, ion, ppm_tolerance):

        # ion backbone score
        if 'ibb' == scoring_alg:
            return ion_backbone_score(spectrum, refseq, ion, ppm_tolerance)

        # ion counting score
        if 'ion' == scoring_alg:
            retindex = 0 if ion == 'b' else 1
            return score_subsequence(spectrum.spectrum, refseq, ppm_tolerance=ppm_tolerance)[retindex]
        
        if 'xcorr' == scoring_alg:
            refspec = gen_spectrum(refseq, ion=ion)['spectrum']
            sparse_ref_spec = make_sparse_array(refspec, .02)
            return xcorr(sparse_spectrum, sparse_ref_spec)

        # default to backbone score
        return backbone_score(spectrum, refseq, ppm_tolerance)

    # hash by the base kmer
    for hittype, hitlist in hits._asdict().items():
        if 'b' in hittype:

            new_binned_b = hash_base_kmers(
                spectrum, 
                hitlist, 
                base_kmer_length, 
                'b'
            )
            
            # add these to the running values
            for k, v in new_binned_b.items():
                base_mer_hashed_b[k] += v

        else:
            new_binned_y = hash_base_kmers(
                spectrum, 
                hitlist, 
                base_kmer_length, 
                'y'
            )
            
            # add these to the running values
            for k, v in new_binned_y.items():
                base_mer_hashed_y[k] += v

    # get rid of any bins that only have 1 hit
    base_mer_hashed_b = {mer: list(set(values)) \
        for mer, values in base_mer_hashed_b.items() \
            if len(list(set(values))) > 1}
    base_mer_hashed_y = {mer: list(set(values)) \
        for mer, values in base_mer_hashed_y.items() \
            if len(list(set(values))) > 1}

    # all of the sequences from b and y hashed kmers that pass filters
    b_seqs, y_seqs = [], []

    # get the overlap score of the sequences from each base kmer
    for _, kmers_b in base_mer_hashed_b.items():
        
        # make it a set
        l = list(set(kmers_b))

        ovscore = overlap_score(l, 'b')

        # if the score is non zero, keep those sequences
        if overlap_score(l, 'b') > 0:

            # score each one and add it to the list
            for kmer_b in l:
                s = score_alg(spectrum, kmer_b, 'b', ppm_tolerance)
                b_seqs = insort_by_index((kmer_b, s), b_seqs, 1)

    # get the overlap score of the sequences from each base kmer
    for _, kmers_y in base_mer_hashed_y.items():

        # make it a set
        l = list(set(kmers_y))

        ovscore = overlap_score(l, 'y')

        # if the score is non zero, keep those sequences
        if overlap_score(l, 'y') > 0:

            # score each one and add it to the list
            for kmer_y in l:
                s = score_alg(spectrum, kmer_y, 'y', ppm_tolerance)
                y_seqs = insort_by_index((kmer_y, s), y_seqs, 1)

    # sort by score
    b_seqs.sort(key=lambda x: x[1], reverse=True)
    y_seqs.sort(key=lambda x: x[1], reverse=True)

    # filter out the good ones
    b_results = [x for x in mean_filtering(b_seqs, mean_filter=2, key=1)]
    y_results = [x for x in mean_filtering(y_seqs, mean_filter=2, key=1)]

    return ([x for x, _ in b_results], [x for x, _ in y_results])