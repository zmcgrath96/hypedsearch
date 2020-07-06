def tryptic(protein: str, missed_cleavages: int) -> list:
    '''
    Take a protien sequence and the number of allowed missed cleavages 
    and return all possible digests of the protein by the missed_cleavages

    Inputs:
        protein:            (str) the protein sequence itself
        missed_cleavages:   (int) the number of missed cleavage sites allowed
    Ouputs:
        (list) all possible digests of the protein at the missed_cleavages
    '''
    # get all the KR indices
    KR_indices = [i for i, aa in enumerate(protein) if aa == 'K' or aa == 'R']

    # if the number of KR indices is less than missed cleavages, then return the protein
    if len(KR_indices) <= missed_cleavages:
        return [protein]

    cleavages = []
    # otherwise go through and take the protein from [startindex+1:end_index+1] to avoid the KR at the start
    # but include it at the end
    for i, idx in enumerate(KR_indices[:-1]):
        end_idx = KR_indices[i + 1 + missed_cleavages] if i + 1 + missed_cleavages < len(KR_indices) else KR_indices[-1]
        cleavages.append(protein[idx+1:end_idx+1])
        
    # take the beginning and end of the protein as well
    start_cleavage = protein[0:KR_indices[missed_cleavages]]
    end_cleavage = protein[KR_indices[-missed_cleavages - 1]:]
    
    return cleavages + [start_cleavage, end_cleavage]