import pandas as pd
from utils.utils import is_gzipped, gunzip_file
###############################################
#               CONSTANTS 
###############################################
FILE_NAME_COL = 'file'
SCORE_COL = 'score'
SCAN_NO_COL = 'scan_no'

score_funcs = ['custom', 'crux']
column_names = {
    'custom': {
        'file_name': 'file',
        'score': 'score',
        'scan_number': 'scan_no'
    },
    'crux': {
        'file_name': 'file',
        'score': 'xcorr score',
        'scan_number': 'scan'
    }
}
###############################################
#             END CONSTANTS
###############################################

###############################################
#           PRIVATE FUNCTIONS
###############################################
def __get_correct_col_names(col_names):
    d_to_l = lambda d: [x for _, x in d.items()]
    for s_func, d in column_names.items():
        ks = d_to_l(d)
        if set.issubset(set(ks), set(col_names)):
            return s_func
    return 'custom'

###############################################
#           END PRIVATE FUNCTIONS
###############################################


'''align_scan_pos

DESC:
    Align scores list to the correspondeing position from the scan numbers
Inputs:
    scores: list of floats of scores
    scan_nos: list of ints of positions where the scores should go
Outputs:
    list, list: aligned scores, aligned scan positions
'''
def align_scan_pos(scores, scan_nos):
    aligned_scores = [0 for _ in range(max(scan_nos) + 1)]
    aligned_scan_nos = [i for i in range(max(scan_nos) + 1)]
    for i, insert_index in enumerate(scan_nos):
        # NOTE: for the mzML files, the first scan starts at 1, not at 0 so minus 1
        # NOTE 2: when i updated to write own scoring the number was off again so removing th e-1
        aligned_scores[insert_index] = scores[i]

    return aligned_scores, aligned_scan_nos

'''get_scores_scan_pos_label

DESC:
    Extract the scores and the scan positions from the file
Inputs:
    file: a string for the filepath to extract results from
kwargs:
    search_substring: a string to search through filenames to limit the search. Defaults to empty
RETUNS:
    list, list: lists of the scores, scan numbers
'''
def get_scores_scan_pos_label(file, search_substring=''):
    sep = '\t' if '.tsv' in file else ','
    if is_gzipped(file):
        file = gunzip_file(file)

    df = pd.read_csv(file, sep, header=0)
        
    s_func = __get_correct_col_names(df.columns)
    col_names = column_names[s_func]

    df = df.sort_values(col_names['score'], ascending=False)
    df = df.drop_duplicates(subset=col_names['scan_number'])
    df = df.sort_values(col_names['scan_number'])
    
    aligned_scores, _ = align_scan_pos(list(df[col_names['score']]), list(df[col_names['scan_number']]))
    return aligned_scores, [], ''

def get_b_y_scores(file: str) -> (list, list):
    '''
    Return both the b and y ion scores

    Inputs:
        file:   string full path to the file name to extract scores from 
    Ouptus:
        (b_scores, y_scores)
        scores both have a list of floats of the scores from each position in order
    '''
    sep = '\t' if '.tsv' in file else ','
    if is_gzipped(file):
        file = gunzip_file(file)

    df = pd.read_csv(file, sep, header=0)
    
    # get the b scores first
    df.sort_values('score_b', ascending=False)
    df.drop_duplicates(subset='scan_no_b')
    df = df.sort_values('scan_no_b')
    b_aligned_scores, _ = align_scan_pos(list(df['score_b']), list(df['scan_no_b']))

    # get the y scores second
    df.sort_values('score_y', ascending=False)
    df.drop_duplicates(subset='scan_no_y')
    df = df.sort_values('scan_no_y')
    y_aligned_scores, _ = align_scan_pos(list(df['score_y']), list(df['scan_no_y']))
    return b_aligned_scores, y_aligned_scores

def pad_scores(score_l1, score_l2, padding=0, side='right'):
    '''
    Makes two lists the same length (filled with 0s)

    Inputs:
        score_l1:   list of floats for the first set of scores
        score_l2:   list of floats for the second set of scores
    kwargs:    
        padding:    float number to pad the list with. Default=0
        side:       string determine which side to pad. Options are {'right', 'left', 'r', 'l'}. Default='right'
    Outputs:
        list, list of modified score_l1, modified score_l2
    '''
    score_l1 = list(score_l1)
    score_l2 = list(score_l2)
    diff = len(score_l1) - len(score_l2)
    if diff == 0:
        return score_l1, score_l2
    elif diff > 0:
        if 'l' in side:
            return score_l1, [padding for _ in range(abs(diff))] + score_l2 
        else: 
            return score_l1, score_l2 + [padding for _ in range(abs(diff))]
    else:
        if 'l' in side: 
            return [padding for _ in range(abs(diff))] + score_l1, score_l2
        else: 
            return score_l1 + [padding for _ in range(abs(diff))], score_l2