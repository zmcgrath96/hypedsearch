from file_io import fasta, csv
from alignment import scoring, filtering

################### Constants ###################
BASE_K = 3
SDEVS = 2
STALL_LENGTH = 3
#################################################

def merge_and_sort(old_b: dict, old_y: dict, new_b: dict, new_y: dict, n=3) -> (dict, dict):
    '''
    Find the best results for each ion type

    Inputs:
        (old_b, old_y, new_b, new_y):   dict where the key indicates the ranking (0 the best) of each subset dictionary.
                                        each dict should have (at least) the entries
                                        {

                                            b_score: float, 
                                            y_score: float, 
                                            sequence: str
                                        }
    kwargs:
        n:  top n number of entries to keep. Default=3
    Outputs:
        (updated_b, updated_y):         (dict, dict) of the same form as the input, but entries sorted
    '''
    to_sort_b = [entry for _, entry in old_b.items()] + [entry for _, entry in new_b.items()]
    to_sort_y = [entry for _, entry in old_y.items()] + [entry for _, entry in new_y.items()]

    # comparing functions for each ion type
    b_cmp = lambda a: a['b_score']
    y_cmp = lambda a: a['y_score']
    sorted_b = sorted(to_sort_b, key=b_cmp, reverse=True)[:n]
    sorted_y = sorted(to_sort_y, key=y_cmp, reverse=True)[:n]

    updated_b, updated_y = {}, {}
    for i in range(n):
        updated_b[i] = sorted_b[i]
        updated_y[i] = sorted_y[i]

    return (updated_b, updated_y)

def search_protein(spectrum: dict, protein_entry: dict, n=3) -> (dict, dict):
    '''
    Search through the protein to find the top n sequences that describe the spectrum

    Inputs:
        spectrum:   dict of the following form:
            {
                scan_no:    int scan number
                spectrum:   list of floats
                level:      int ms level
            }
        database:   dict from a .fasta file with the form
            {
                sequence: str,
                name: str, 
                id: str
            }
    kwargs:
        n:      int top number of results to return for both b and y ion scores
    Outptus:
        (b_dict, y_dict)
        Both of these ion dictionaries have the top n entries keyed by their rank (0 to n-1, 0 being the best)
         {
             0: {starting_position: int, ending_position: int, k: int, b_score: float, y_score: float, sequence: str, protein_name: str, protein_id: str}
             ...
             n-1: {...}
         }
    '''
    # get any significant kmers from the scoring
    base_scores = scoring.score_sequence(spectrum['spectrum'], protein_entry['sequence'], BASE_K)
    b_anchors = filtering.stddev_filter(base_scores, 'b_score', SDEVS)
    y_anchors = filtering.stddev_filter(base_scores, 'y_score', SDEVS)
    # extend these interesting kmers as much as possible
    (extended_bs, extended_ys) = scoring.kmer_extend(spectrum['spectrum'], protein_entry['sequence'], b_anchors, y_anchors, STALL_LENGTH)
    # take the top n best kmers
    b_res, y_res = {}, {}
    for i in range(n):
        b_res[i] = extended_bs[i]
        b_res[i]['protein_name'] = protein_entry['name']
        b_res[i]['protein_id'] = protein_entry['id']
        y_res[i] = extended_ys[i]
        y_res[i]['protein_name'] = protein_entry['name']
        y_res[i]['protein_id'] = protein_entry['id']

    return b_res, y_res

def search_proteins(spectrum: dict, database: list, n=3) -> (dict, dict):
    '''
    Search through the proteins to find the top n sequences that describe the spectrum

    Inputs:
        spectrum:   dict of the following form:
            {
                scan_no:    int scan number
                spectrum:   list of floats
                level:      int ms level
            }
        database:   list of dicts from a .fasta file with the form
            {
                sequence: str,
                name: str, 
                id: str
            }
    kwargs:
        n:      int top number of results to return for both b and y ion scores
    Outptus:
        (b_dict, y_dict)
        Both of these ion dictionaries have the top n entries keyed by their rank (0 to n-1, 0 being the best)
         {
             0: {starting_position: int, ending_position: int, k: int, b_score: float, y_score: float, sequence: str, protein_name: str}
             ...
             n-1: {...}
         }
    '''
    b_results, y_results = {}, {}
    # go through every protein in the database
    for protein_entry in database:
        top_prot_b, top_prot_y = search_protein(spectrum, protein_entry, n)
        # if b_results and y_results are empty, just set them to the results we have for the first protein
        if top_prot_b is None or len(b_results.keys()) == 0:
            b_results = top_prot_b
        if top_prot_y is None or len(y_results.keys()) == 0:
            y_results = top_prot_y

        b_results, y_results = merge_and_sort(b_results, y_results, top_prot_b, top_prot_y, n)

    return (b_results, y_results)
