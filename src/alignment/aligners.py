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

    b_prots = [x['protein_name'] for _, x in b_scores.items()]
    y_prots = [x['protein_name'] for _, x in y_scores.items()]
    # see if any pairing does exist
    if not (any([True if x in y_prots else False for x in b_prots])):
        return pairings

    # at least 1 pairing exists
    for rank, bscore in b_scores.items():
        prot = bscore['protein_name']
        # if the protein is not found in y proteins, skip
        if prot not in y_prots:
            continue
        # if the protein is not in pairings, add it
        if prot not in pairings: 
            pairings[prot] = []
        # find the pairing we know exists and add it
        for yrank, yscore in y_scores.items():
            if yscore['protein_name'] == prot:
                pairings[prot].append((rank, yrank, bscore, yscore))

    return pairings


def align_spectrum_by_protein_ions(spectrum: list, b_scores: dict, y_scores: dict, n=3) -> dict:
    '''
    CREATED 15 APRIL 2020
    Create an alignment from b and y scores. If any alignments can be made, then they will be from the same protein

    Inputs:
        spectrum:               list of floats the mass spectrum to which the scores are to be aligned
        (b_score, y_score):     dict of top ion scores with the form 
            {
                0: {starting_position: int, ending_position: int, b_score: float, y_score: float, sequence: str, protein_name: str}
                ...
                n-1: {...}
            }
    kwargs:
        n:          int number of alignments to return. If n results aren't available, all results returned. Default=3
    Outputs:
        dict        dictionary of results with the key being the index of the ranking (0 being the highest rank) of the form
            {
                0: {starting_position: int, ending_position: int, length: int, b_score: float, y_score: float, confidence: float, protein_name: string, spectrum: list},
                ...
                n-1: {...}
            }
    '''
    pairings = find_protein_pairings(b_scores, y_scores)
    pairings_with_confidence = []
    for protname, pairing_list in pairings.items():
        for pairing in pairing_list:
            pairings_with_confidence.append((scoring.confidence_simple(pairing[2], pairing[3]), protname, pairing[2], pairing[3]))

    pairings_with_confidence.sort(key=lambda x: x[0], reverse=True)
    alignments = {}
    num_alignments = min(n, len(pairings_with_confidence))
    for i in range(num_alignments):
        p = pairings_with_confidence[i]
        alignments[i] = {
            'starting_position': p[2]['starting_position'],
            'ending_position': p[3]['ending_position'],
            'length': p[3]['ending_position'] - p[2]['starting_position'] + 1,
            'b_score': p[2]['b_score'],
            'y_score': p[3]['y_score'],
            'confidence': p[0], 
            'protein_name': p[1],
            'spectrum': spectrum
        }

    return alignments

    