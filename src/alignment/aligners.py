from src.alignment import scoring

def find_protein_pairings(b_scores: dict, y_scores: dict) -> dict:
    '''
    Find any pairings of the b and y score by parent proteins.

    Inputs:
        (b_score, y_score):    dict of top ion scores with the form 
            {
                0: {starting_position: int, ending_position: int, k: int, b_score: float, y_score: float, sequence: str, protein_name: str}
                ...
                n-1: {...}
            }
    Outputs:
        dict        dictionary of the pairings of the form
            {
                protein_name: [
                    (b_rank, y_rank, b_score, y_score)  
                ]
            }
    '''
    pairings = {}

    b_prots = [x['protein_name'] for _, x in b_scores]
    y_prots = [x['protein_name'] for _, x in y_scores]
    # see if any pairing does exist
    if not (any([True if x in y_prots else False for x in b_prots])):
        return pairings

    # at least 1 pairing exists
    for rank, bscore in b_scores:
        prot = bscore['protein_name']
        # if the protein is not found in y proteins, skip
        if prot not in y_prots:
            continue
        # if the protein is not in pairings, add it
        if prot not in pairings: 
            pairings[prot] = []
        # find the pairing we know exists and add it
        for yrank, yscore in y_scores:
            if yscore['protein_name'] == prot:
                pairings[prot].append((rank, yrank, bscore, yscore))

    return pairings


def align_ion_scores_same_prot(b_score: dict, y_score: dict, position_thresh=None) -> dict:
    '''
    Try and find an alignment between 2 ion scores. These ions should have the same protein

    Inputs:
        (b_score, y_score): dictionaries with the entries
            {starting_position: int, ending_position: int, k: int, b_score: float, y_score: float, sequence: str, protein_name: str}
    kwargs: 
        position_thresh:    int number of positions each ion alignment is allowed to be from the other before 
                            no attempt at alignment is made. If none, an alignment is attempted for the entire protein. Defautl=None
    Outputs:
        dict        dictionary with information on an alignemnt. If none is found, then empty dict returned
            {'starting_position': int, 'ending_position': int, 'b_score': float, 'y_score': float, 'confidence': float}                    
    '''


def align_spectrum_from_ion_scores(spectrum: list, b_scores: dict, y_scores: dict, n=3) -> dict:
    '''
    CREATED 15 APRIL 2020
    Create an alignment from b and y scores. Total score for now is the sum of the two scores that best identify the sequence

    Inputs:
        spectrum:               list of floats the mass spectrum to which the scores are to be aligned
        (b_score, y_score):     dict of top ion scores with the form 
            {
                0: {starting_position: int, ending_position: int, k: int, b_score: float, y_score: float, sequence: str, protein_name: str}
                ...
                n-1: {...}
            }
    kwargs:
        n:          int number of alignments to return. If n results aren't available, all results returned. Default=3
    Outputs:
        dict        dictionary of results with the key being the index of the ranking (0 being the highest rank) of the form
            {
                0: {starting_position: int, ending_position: int, length: int, score: float, protein: string},
                ...
                n-1: {...}
            }
    '''
    pairings = find_protein_pairings(b_scores, y_scores)

    