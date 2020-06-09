from src.types.objects import KmerMassesResults, Spectrum
from src.scoring.scoring import score_subsequence
from src.utils import insort_by_index

from statistics import mean
from typing import Iterable
from collections import defaultdict

def slope_filtering(a: Iterable, min_window_size=5, mean_filter=1) -> list:
    '''
    Filter out values by slope. Use a window size of at least min_window_size and at max 1% the size of 
    iterable a. Result is the subset of results whose window mean > mean_filter*mean_slope. a should be sorted
    with the steepest slope at the earliest indices 
    
    Inputs:
        a:     (Iterable) the list or list like object to filter
    kwargs:
        min_window_size:    (int) the minimum window size to use for sliding means. Default=5
        mean_filter:        (int) a scaling factor to determine what values pass. This value
                                  is mutliplied by the mean slope and values > than this are accepted. Default=1
    Outputs:
        list with the filtered results
    '''
    if len(a) < min_window_size:
        return a
    adjusted_window_size = len(a) // 100
    window_size = min_window_size if adjusted_window_size < min_window_size else adjusted_window_size
    slopes = [(a[i+window_size] - a[i])/window_size for i in range(len(a) - window_size)]
    avgslope = mean(slopes)
    filtered = []
    for i in range(len(slopes)):
        if abs(slopes[i]) > abs(avgslope):
            filtered.append(a[i])
        else:
            break
    return filtered


def result_filtering(spectrum: Spectrum, hits: KmerMassesResults, base_kmer_length: int):
    # narrow down the results to a few promising kmers
    basemerhashedb = defaultdict(list)
    basemerhashedy = defaultdict(list)
    basemerbblacklist = {}
    basemeryblacklist = {}
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
                    basemerbscore = score_subsequence(spectrum.spectrum, basemerb)[0]
                    # blacklist for bad scores
                    if not basemerbscore > 0:
                        basemerbblacklist[basemerb] = None
                        continue
                    # insert in order
                    b_scores = insort_by_index((basemerb, score_subsequence(spectrum.spectrum, basemerb)[0]), b_scores, 1)
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
                    basemeryscore = score_subsequence(spectrum.spectrum, basemery)[1]
                    # blacklist for bad score
                    if not basemeryscore > 0:
                        basemeryblacklist[basemery] = None
                        continue
                    # insert in order
                    y_scores = insort_by_index((basemery, score_subsequence(spectrum.spectrum, basemery)[1]), y_scores, 1)
                basemerhashedy[basemery].append(masssequence_hit)

    toscoreb = []
    toscorey = []
    
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
    
    best_b_results = sorted(toscoreb, key=lambda mer: score_subsequence(spectrum.spectrum, mer)[0], reverse=True)
    best_y_results = sorted(toscorey, key=lambda mer: score_subsequence(spectrum.spectrum, mer)[1], reverse=True)

    return (best_b_results, best_y_results)