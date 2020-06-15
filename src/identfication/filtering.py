from src.types.objects import KmerMassesResults, Spectrum
from src.scoring.scoring import score_subsequence, backbone_score, ion_backbone_score
from src.utils import insort_by_index

from statistics import mean
from typing import Iterable
from collections import defaultdict

def slope_filtering(a: Iterable, min_window_size=5, mean_filter=1) -> list:
    '''
    Filter out values by slope. Use a window size of at least min_window_size and at max 1% the size of 
    iterable a. Result is the subset of results whose window mean > mean_filter*mean_slope. a should be sorted
    with the steepest slope at the earliest indices 

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
    
    # calculate the slope for a window size
    slopes = [(a[i+window_size] - a[i])/window_size for i in range(len(a) - window_size)]
    avgslope = mean(slopes)
    filtered = []
    # keep anything with a slope > filter
    for i in range(len(slopes)):
        if abs(slopes[i]) > abs(avgslope):
            filtered.append(a[i])
        else:
            break
    return filtered


def result_filtering(spectrum: Spectrum, hits: KmerMassesResults, base_kmer_length: int, ppm_tolerance=20, scoring_alg='ibb') -> (list, list):
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
    basemerhashedb = defaultdict(list)
    basemerhashedy = defaultdict(list)
    
    # scoring algorithm to use
    def score_alg(refseq, ion):
        if 'ibb' == scoring_alg:
            return ion_backbone_score(spectrum, refseq, ion, ppm_tolerance)
        elif 'ion' == scoring_alg:
            retindex = 0 if ion == 'b' else 1
            return score_subsequence(spectrum.spectrum, refseq, ppm_tolerance=ppm_tolerance)[retindex]
        else:
            return backbone_score(spectrum, refseq, ppm_tolerance)


    # keep track of base kmers that don't score > 0
    basemerbblacklist = {}
    basemeryblacklist = {}

    # sorted lists of (basekmers, score) by score
    b_scores = []
    y_scores = []

    # hash by the base kmer
    for hittype, hitlist in hits._asdict().items():
        if 'b' in hittype:
            for masssequence_hit in hitlist:
                basemerb = masssequence_hit[:base_kmer_length]
                # check to see if we've seen it
                # skip if its blacklisted
                if basemerb in basemerbblacklist:
                    continue
                # if we've not seen it, try and add it
                if basemerb not in basemerhashedb:
                    basemerbscore = score_alg(basemerb, 'b')
                    # blacklist for bad scores
                    if not basemerbscore > 0:
                        basemerbblacklist[basemerb] = None
                        continue
                    # insert in order
                    b_scores = insort_by_index((basemerb, basemerbscore), b_scores, 1)
                basemerhashedb[basemerb].append(masssequence_hit)
        else:
            for masssequence_hit in hitlist:
                basemery = masssequence_hit[len(masssequence_hit)-base_kmer_length:]
                # check to see if we've seen it
                # skip if its blacklisted
                if basemery in basemeryblacklist:
                    continue
                # if we've not seen it, try and add it
                if basemery not in basemerhashedy:
                    basemeryscore = score_alg(basemery, 'y')
                    # blacklist for bad score
                    if not basemeryscore > 0:
                        basemeryblacklist[basemery] = None
                        continue
                    # insert in order
                    y_scores = insort_by_index((basemery, basemeryscore), y_scores, 1)
                basemerhashedy[basemery].append(masssequence_hit)

    toscoreb = []
    toscorey = []
    
    # mean filtered base kmers by scores
    b_filtered = [b_scores[-1*i - 1] for i in range(len(slope_filtering([x[1] for x in b_scores[::-1]])))]
    y_filtered = [y_scores[-1*i - 1] for i in range(len(slope_filtering([x[1] for x in y_scores[::-1]])))]
    
    if len(b_filtered) < 5: # default to be able to report some result
        b_filtered = b_scores[::-1][:5]
    if len(y_filtered) < 5:
        y_filtered = y_scores[::-1][:5]
    
    # revese to go from highest to lowest
    for basemerb in b_filtered:
        toscoreb += [basemerb[0]] + basemerhashedb[basemerb[0]]
        
    for basemery in y_filtered:
        toscorey += [basemery[0]] + basemerhashedy[basemery[0]]
    
    # sorted score of all leftover sequences from highest to lowest
    best_b_results = sorted(toscoreb, key=lambda mer: score_alg(mer, 'b'), reverse=True)
    best_y_results = sorted(toscorey, key=lambda mer: score_alg(mer, 'y'), reverse=True)

    return (best_b_results, best_y_results)