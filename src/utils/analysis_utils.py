from numpy import argmax, average, argsort
from math import inf
from scipy.signal import find_peaks
from typing import List, Dict

'''__get_argmax_max

DESC:
    find the argmax and the max of an iterable
Inputs:
    i: iterable to return max and argmax on
Outputs:
    int, val: argmax, max
'''
def __get_argmax_max(i):
    return argmax(i), max(i)

'''get_highest_scoring

DESC:
    return the list and key of the best scoring list of scores
Inputs:
    scores: iterable of lists each entry of the lists should be a number. If scores is not an iterable, [], 0 is returned
kwargs:
    measure: string determines which metric is used for the calculate best score. Can be 'max', 'average', 'sum'. Default=average
Outputs:
    list, val: the list and the index or key of the best score list
'''
def get_highest_scoring(scores, measure='average'):
    func = sum if 'sum' in measure.lower() else (max if 'max' in measure.lower() else average)
    
    l = None
    v = None
    if type(scores) is dict:
        for key, score in scores.items():
            if l is None:
                l = score
                v = key
            if func(score) > func(l):
                l = score
                v = key
    else:
        for i, score in enumerate(scores):
            if l is None:
                l = score
                v = i
            if func(score) > func(l):
                l = score 
                v = i
    return l, v

'''get_top_n

DESC:
    returns the top n list and keys or indices of an interable of lists
Inputs:
    scores: iterable of lists each entry of the lists should be a number. If scores is not an iterable, [] is returned
kwargs:
    measure: string determines which metric is used for the calculate best score. Can be 'max', 'average', 'sum'. Default=average
    n: int number of top scores to return. Default=5
Outputs:
    [(list, val)]: list of tuples with the list, index or key of the top n in order
'''
def get_top_n(scores, measure='average', n=5):
    top_n = []
    n = len(scores) if len(scores) < n else n

    for _ in range(n):
        l, v = get_highest_scoring(scores, measure=measure)
        t = (l, v)
        top_n.append(t)
        del scores[v]
    
    return top_n

'''__find_top_n_peaks

DESC: 
    finds the top n peaks of a signal 
Inputs:
    signal: a list of floats 
kwargs: 
    n: top number of peaks to find. Default=5
    height: float lowest number to consider for peaks
Outputs:
    list of ints. These are the positions of the top peaks
'''
def __find_top_n_peaks(signal, n=5, height=0):
    peaks, _ = find_peaks(signal, height=height)
    # returns indices that would sort it lowest to hightest
    peak_w_values = [(x, signal[x]) for x in peaks]
    peak_w_values.sort(key = lambda x: x[1])
    top_n = peak_w_values[-n:]
    top_n_indx = [x[0] for x in top_n]
    return top_n_indx 

def get_top_n_prots(prots: dict, n=5) ->List[Dict]:
    '''
    find the top n scoring proteins

    Inputs:
        prots: dictionary of scores where the keys are protein names
    kwargs: 
        n: int number of top scores to find. Default=5
    Outputs:
        List of dictionaries [{protein_name:str, score: float, position: int}]
    '''
    if type(prots) is not dict:
        raise Exception('prots should be a dictionary. Type is {}'.format(type(prots)))
    
    top_n_each = []
    for prot_name, score in prots.items():
        top_n = __find_top_n_peaks(score)
        top_n_each += [{'protein_name': prot_name, 'score': score[i], 'position': int(i)} for i in top_n]

    top_n_each.sort(key = lambda x: x['score'], reverse=True)
    t = top_n_each[:n]
    return t