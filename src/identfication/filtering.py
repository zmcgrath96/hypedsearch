from src.types.objects import KmerMassesResults, Spectrum
from src.scoring.scoring import score_subsequence, backbone_score, ion_backbone_score, intensity_ion_backbone_score, xcorr
from src.utils import insort_by_index, make_sparse_array
from src.spectra.gen_spectra import gen_spectrum

from statistics import mean, stdev
from typing import Iterable
from collections import defaultdict
from operator import itemgetter

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

def stddev_filtering(a: Iterable, stddevs=1, key=None) -> list:
    '''
    Filter out values that fall < stddevs * stddev. Key allows usage of stuctures in the iterable.

    Inputs:
        a:          (Iterable) values to sort
    kwargs:
        stddevs:    (float) the number of standard deviations to use as the filter. Default=1
        key:        (any) the key to use to access values in elements. If None, it is assumed that
                            it is an iterable of numbers. Default=None
    Outputs:
        (list) the filtered list of values
    '''
    # function for getting values
    get_val = lambda x: x if key is None else (x[key] if type(x) == dict or type(x) == list or type(x) == tuple else getattr(x, key))

    # get the standard deviation
    stddev = stdev([get_val(e) for e in a])

    # return anything above the filter
    return [e for e in a if get_val(e) > stddevs * stddev]

def hash_base_kmers(
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

def overlap_assembler(l: list, ion: str) -> str:
    '''
    Take a list of sequences and return the longest sequence that contains
    the most of the rest of the sequences.

    Example:
        l: [ABCDEFG, ABCDE, ABCXY, ABCD, ABC]

        ABCDEFG contains [ABCDE, ABCD, ABC] which is the most of any others
        so ABCDEFG would be returned

    Inputs:
        l:      (list) the sequences to assemble together   
        ion:    (str) the ion this is performed for. y goes right to left, b left to right
    Outputs:
        (str) the sequence that contains the most of the others. If no
                sequence can be assembled, None is returned
    '''
    # sort longest to shortest
    l.sort(key=len, reverse=True)
    
    # got through each element and count the number of smaller 
    # sequences it contains
    assembler = defaultdict(lambda: 0)
    for i, e in enumerate(l[:-1]):
        for e2 in l[i+1:]:
            
            # get left to right or right to left depending on the ion
            cmp_str = e[:len(e2)] if ion == 'b' else e[-len(e2):]
                        
            if e2 == cmp_str:
                assembler[e] += 1
                
    # return the str if possible, None otherwise
    try:    
        max_count = max([(k, v) for k, v in assembler.items()], key=itemgetter(1))
        return max_count[0] if max_count[1] > 1 else None
    except:
        return None


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

        if 'iibb' == scoring_alg:
            return intensity_ion_backbone_score(spectrum, refseq, ion, ppm_tolerance) 

        # default to backbone score
        return backbone_score(spectrum, refseq, ppm_tolerance)

    # hash by the base kmer
    for hittype, hitlist in hits._asdict().items():
        if 'b' in hittype:

            new_binned_b = hash_base_kmers(
                hitlist, 
                base_kmer_length, 
                'b'
            )
            
            # add these to the running values
            for k, v in new_binned_b.items():
                base_mer_hashed_b[k] += v

        else:
            new_binned_y = hash_base_kmers(
                hitlist, 
                base_kmer_length, 
                'y'
            )
            
            # add these to the running values
            for k, v in new_binned_y.items():
                base_mer_hashed_y[k] += v

    # get rid of any bins that only have 1 hit
    base_mer_hashed_b = {mer: values \
        for mer, values in base_mer_hashed_b.items() \
            if len(list(set(values))) > 1}

    base_mer_hashed_y = {mer: values \
        for mer, values in base_mer_hashed_y.items() \
            if len(list(set(values))) > 1}

    # all of the sequences from b and y hashed kmers that pass filters
    b_seqs, y_seqs = [], []

    def reduce_sequences(d: dict, ion: str):
        seqs = []
        v: list
        for k, v in d.items():
       
            reduced_seq = overlap_assembler(list(set(v)), ion)

            if reduced_seq is None:
                seqs.append(str(k))
                continue

            seqs.append(reduced_seq)
        return seqs

    def reduce_sequences2(d: dict, ion: str):
        seqs = []
        v: list
        for k, v in d.items():
            reduced_seq = max(
                [(seq, score_alg(spectrum, seq, ion, ppm_tolerance)) for seq in list(set([str(k)] + v))],
                key=itemgetter(1)
            )[0]
            seqs.append(reduced_seq)
        return seqs

    # print(f'Base mer hashed b:\n{base_mer_hashed_b}')
    # print(f'Base mer hashed y:\n{base_mer_hashed_y}')

    # reduce to the most overlapping sequences
    b_seqs = reduce_sequences(base_mer_hashed_b, 'b')
    y_seqs = reduce_sequences(base_mer_hashed_y, 'y')

    # score them and take non zero scores
    b_results = [(x, score_alg(spectrum, x, 'b', ppm_tolerance)) \
        for x in b_seqs if score_alg(spectrum, x, 'b', ppm_tolerance) > 0]
    y_results = [(x, score_alg(spectrum, x, 'y', ppm_tolerance)) \
         for x in y_seqs if score_alg(spectrum, x, 'y', ppm_tolerance) > 0]

    # sort the results by score high to low
    b_results.sort(key=itemgetter(1), reverse=True)
    y_results.sort(key=itemgetter(1), reverse=True)

    # print(f'B results before filtering:\n{b_results}')
    # print(f'Y results before filtering:\n{y_results}')

    # take the scores that pass our filter
    def filter_scores(l: list) -> list:
        if len(l) < 2:
            return l

        filtered = [
            stddev_filtering(l, stddevs=2, key=1),
            mean_filtering(l, mean_filter=2, key=1),
            slope_filtering(l, mean_filter=2, key=1)
        ]
        
        nonzero = [x for x in filtered if len(x) > 0]
        if len(nonzero) == 0:
            return []
        
        return min(nonzero, key=len)


    filtered_b_results = filter_scores(b_results)
    filtered_y_results = filter_scores(y_results)

    # if we have nothing, take the top 5 scores
    if len(filtered_b_results) == 0:

        top_scores = []
        
        # b results is sorted, so go through and just get the top 5
        for b_res in b_results:
            
            # get the score and see if its in the list
            if b_res[1] not in top_scores and len(top_scores) < 5:
                top_scores.append(b_res[1])
            
            # if the list is at 5, break
            elif b_res[1] not in top_scores and len(top_scores) >= 5:
                break

            # add it to the list otherwise
            filtered_b_results.append(b_res)

    # if we have nothing, take the top 5 scores
    if len(filtered_y_results) == 0:

        top_scores = []
        
        # b results is sorted, so go through and just get the top 5
        for y_res in y_results:
            
            # get the score and see if its in the list
            if y_res[1] not in top_scores and len(top_scores) < 5:
                top_scores.append(y_res[1])
            
            # if the list is at 5, break
            elif y_res[1] not in top_scores and len(top_scores) >= 5:
                break

            # add it to the list otherwise
            filtered_y_results.append(y_res)

    return ([x for x, _ in filtered_b_results], [x for x, _ in filtered_y_results])