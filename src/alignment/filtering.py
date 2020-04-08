import numpy as np

def stddev_filter(scores: list, score_key: str, stddevs=2) -> list:
    '''
    Filter scores by a number of standard deviations

    Inputs:
        scores:     list of dictionaries
        score_key:  str key value to index into scores to get the numeric score value
    kwargs:
        stddevs:    int number of standard deviations to use as the cuttoff 
    Outputs:
        list        a subset of the scores list with the filtered scores        
    '''
    stddev = np.std([score[score_key] for score in scores])
    filtered = [score for score in scores if score[score_key] >= stddev*stddevs]
    return filtered