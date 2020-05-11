import numpy as np

def stddev_filter(scores: list, score_key: str, stddevs=2) -> list:
    '''
    Filter scores by a number of standard deviations

    Inputs:
        scores:     list of dictionaries or objects
        score_key:  (str) key value to index into scores to get the numeric score value
    kwargs:
        stddevs:    (int) number of standard deviations to use as the cuttoff 
    Outputs:
        list        a subset of the scores list with the filtered scores        
    '''
    getval = lambda x: getattr(x, score_key) if type(scores[0]) != dict else lambda x: x[score_key]
    stddev = np.std([getval(score) for score in scores])
    filtered = [score for score in scores if getval(score) >= stddev*stddevs]
    return filtered