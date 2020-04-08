from alignment import comparisons
from spectra import gen_spectra

def __new_entry(old_entry: dict, prot: str, spectrum: list, ion='b') -> dict:
    '''
    Generate a new entry from the old entry
    
    Input:
        old_entry:   dict entry with k, sequence, b and y scores, start and end positions
        prot:        str sequence of the protein
        spectrum:    list spectrum to score against
    kwargs:
        ion:         str ion type to determine which. Options are 'b', 'y'. Default='b' 
    Ouptut:
        new_entry:   dict entry with the new k, new sequence, new b and y scores, new start and end positions
    '''
    starting_pos = old_entry['starting_position'] if ion == 'b' else old_entry['starting_position'] - 1
    ending_pos = old_entry['ending_position'] + 1 if ion == 'b' else old_entry['ending_position']
    if starting_pos < 0 or ending_pos > len(prot) - 1:
        return old_entry
    mer_seq = prot[starting_pos:ending_pos+1]
    mer_spec_b = gen_spectra.gen_spectrum(mer_seq, ion='b')['spectrum']
    mer_spec_y = gen_spectra.gen_spectrum(mer_seq, ion='y')['spectrum']
    return {
        'k': old_entry['k'] + 1,
        'sequence': mer_seq,
        'starting_position': starting_pos,
        'ending_position': ending_pos,
        'b_score': comparisons.compare_masses(spectrum, mer_spec_b),
        'y_score': comparisons.compare_masses(spectrum, mer_spec_y)
    }

def kmer_extend(spectrum: list, sequence: str, b_kmers: list, y_kmers: list, stall_length=3) -> (list, list):
    '''
    Extend kmers until the score does not increase any more. Stall length is used to determine
    the number of iterations the score can go without increasing without finishing the kmer growth

    Inputs:
        spectrum:           list of floats the mass spectrum being checked
        sequence:           str the sequence to add amino acids from
        (b_kmers, y_kmers): list of dictionaries with (at least) the following indices
                        {
                            b_score: float,
                            y_score: float
                            k: int,
                            starting_position: int,
                            ending_position: int
                        }
    kwargs:
        stall_length:   int the number of iterations a subsequence is allowed to go witth 
                        no increase in score before finishing kmer growth on a certain kmer
    Outputs:
        (b_list, y_list)
        list of reversed order rankings of best matches for each ion type. Index 0 is the best match and index (len-1) is the worst
    '''
    # keep track of the ones that we want to increment. Prepend it to the front and take the top
    b_list = []
    y_list = []
    b_stall_counter = {b['starting_position']: [stall_length, b] for b in b_kmers}
    y_stall_counter = {y['ending_position']: [stall_length, y] for y in y_kmers}

    # go through the kmers and extend them for as long as possible
    while(len(b_kmers)):
        b_tmp = []
        # rounds of incrementing
        for i in range(len(b_kmers)):
            # try extending the current entry
            updated = __new_entry(b_kmers[i], sequence, spectrum, ion='b')
            # if the score increased, keep it
            if updated['b_score'] > b_kmers[i]['b_score']:
                b_tmp.append(updated)
                b_stall_counter[b_kmers[i]['starting_position']][1] = updated
            # if the score didnt increase and the counter is non zero, decrement the counter
            elif b_stall_counter[b_kmers[i]['starting_position']][0] > 0:
                b_tmp.append(updated)
                b_stall_counter[b_kmers[i]['starting_position']][0] -= 1
            # the entry has run out of chances
            else:
                b_list.insert(0, b_kmers[i])
        b_kmers = b_tmp

    # go through the kmers and extend them for as long as possible  
    while(len(y_kmers)):
        y_tmp = []
        # rounds of incrementing
        for i in range(len(y_kmers)):
            # try extending the current entry
            updated = __new_entry(y_kmers[i], sequence, spectrum, ion='y')
            # if the score increased, keep it
            if updated['y_score'] > y_kmers[i]['y_score']:
                y_tmp.append(updated)
                y_stall_counter[y_kmers[i]['ending_position']][1] = updated
            # if the score didnt increase and the counter is non zero, decrement the counter
            elif y_stall_counter[y_kmers[i]['ending_position']][0] > 0:
                y_tmp.append(updated)
                y_stall_counter[y_kmers[i]['ending_position']][0] -= 1
            # the entry has run out of chances
            else: 
                y_list.insert(0, y_kmers[i])
        y_kmers = y_tmp
        
    return (b_list, y_list)

def score_subsequence(pepspec: list, subseq: str) -> (float, float):
    '''
    Score a mass spectrum to a substring of tagged amino acids

    Inputs:
        pepspec:    list of floats the mass spectrum to score
        subseq:     str of tagged amino acids to score the spectrum against
    Outputs:
        (b_score, y_score): (float, float) the b and y ion scores generated from this comparison
    '''
    kmerspec_b = gen_spectra.gen_spectrum(subseq, ion='b')['spectrum']
    kmerspec_y = gen_spectra.gen_spectrum(subseq, ion='y')['spectrum']
    b_score = comparisons.compare_masses(pepspec, kmerspec_b)
    y_score = comparisons.compare_masses(pepspec, kmerspec_y)
    return (b_score, y_score)

def score_sequence(pepspec: list, prot: str, k: int) -> (list, list):
    '''
    Generate scores of a spectrum against a protein kmers

    Inputs:
        pepspec:    list of floats the mass spectrum in question
        prot:       str the string sequence of amino acids of a protein
        k:          int the kmer size to break the protein into
    Outputs:
        list of dictionaries of the form
        {
            'k': int,
            'sequence': str, 
            'b_score': float,
            'y_score': float,
            'starting_position': int,
            'ending_position': int
        }
    '''
    kmer_scores = []
    make_mer_sp = lambda starting_pos, prot, k: prot[starting_pos: starting_pos+k]
    for i in range(len(prot) - k + 1):
        kmer = make_mer_sp(i, prot, k)
        b_score, y_score = score_subsequence(pepspec, kmer)
        entry = {
            'k': k,
            'sequence': kmer,
            'b_score': b_score,
            'y_score': y_score,
            'starting_position': i, 
            'ending_position': i + k - 1
        }
        kmer_scores.append(entry)

    return kmer_scores