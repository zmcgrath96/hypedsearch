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

def score_filter(scores: list, score_key: str, score=0, metric='g') -> list:
    '''
    Filter scores by value

    Inputs:
        scores:     (list) dictionaries or object with score keys
        score_key:  (str) key with which to index the scores
    kwargs:
        score:      (float) the value to use as the filter. Default=0
        metric:     (str) Use Greater than (g), Greater than or equal to (ge), 
                    Less than (l), Less than or equal to (le), or Equal to (e). Possible values [g, ge, l, le, e]. Default='g'
    Outputs:
        (list) subset of the scores list passed in
    '''
    comp = lambda a, b: a > b if metric.lower() == 'g' else (lambda a, b: a >= b if metric.lower() == 'ge' else \
            (lambda a, b: a == b if metric.lower() == 'e' else (lambda a, b: a <= b if metric.lower() == 'le' else \
                lambda a, b: a < b)))
    getscore = lambda x: x[score_key] if type(scores[0]) == dict else getattr(x, score_key)
    return [s for s in scores if comp(getscore(s), score)]