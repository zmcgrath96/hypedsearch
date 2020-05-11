from src.scoring import mass_comparisons
from src.spectra import gen_spectra

def new_entry(old_entry: dict, prot: str, spectrum: list, ion='b') -> dict:
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
    # check that we are operating on a valid entry
    keys = ['k', 'starting_position', 'ending_position', 'b_score', 'y_score']
    if any([old_entry[k] is None for k in keys]) or prot is None or spectrum is None:
        return old_entry

    starting_pos = old_entry['starting_position'] if ion == 'b' else old_entry['starting_position'] - 1
    ending_pos = old_entry['ending_position'] + 1 if ion == 'b' else old_entry['ending_position']
    if starting_pos < 0 or ending_pos > len(prot) - 1:
        return old_entry

    # check for negative lengths
    if starting_pos > ending_pos or ending_pos < starting_pos:
        return old_entry

    mer_seq = prot[starting_pos:ending_pos+1]
    mer_spec_b = gen_spectra.gen_spectrum(mer_seq, ion='b')['spectrum']
    mer_spec_y = gen_spectra.gen_spectrum(mer_seq, ion='y')['spectrum']
    return {
        'k': old_entry['k'] + 1,
        'sequence': mer_seq,
        'starting_position': starting_pos,
        'ending_position': ending_pos,
        'b_score': mass_comparisons.compare_masses(spectrum, mer_spec_b),
        'y_score': mass_comparisons.compare_masses(spectrum, mer_spec_y)
    }

def extend_kmer(spectrum: list, sequence: str, kmer: dict, ion: str, stall_length=3) -> dict:
    '''
    Extend a kmer until the score tells us that the adding amino acids doens't make it a better alignment
    
    Inputs:
        spectrum:       list of floats. The mass spectrum in question
        sequence:       str The full protein sequence we are pulling amino acids from 
        kmer:           dict of the form
                        {
                            b_score: float, 
                            y_score: float,
                            k: int, 
                            starting_position: int,
                            ending_position: int,
                        }
        ion:            str the ion type we are looking at. Should be 'b' or 'y'
    kwargs:
        stall_length:   int the number of iterations a subsequence is allowed to go witth 
                        no increase in score before finishing kmer growth on a certain kmer. Default=3
    Outputs
        dict with updated values of the form
            {
                b_score: float,
                y_score: float, 
                k: int, 
                starting_position: int,
                ending_positoin: int,
            }
    '''
    if ion.lower() not in ['b', 'y']:
        return kmer
    score_key = 'b_score' if ion.lower() == 'b' else 'y_score'
    # keep track of the last time a score increased
    last_maintenance = kmer
    # keep going until we run out of extension
    while stall_length > 0:
        updated = new_entry(kmer, sequence, spectrum, ion=ion)
        if updated[score_key] > kmer[score_key] and updated[score_key] > 0 and updated['k'] != kmer['k']:
            last_maintenance = updated
        else: 
            stall_length -= 1
        kmer = updated
    return last_maintenance