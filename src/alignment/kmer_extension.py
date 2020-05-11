from src.scoring import mass_comparisons
from src.spectra import gen_spectra
from src.modules.kmer import Kmer 
from src.modules.spectrum import Spectrum
from src.modules.scores import Scores

def new_entry(old_entry: Scores, prot: str, spectrum: Spectrum, ion='b') -> Scores:
    '''
    Generate a new entry from the old entry
    
    Input:
        old_entry:   Scores object
        prot:        str sequence of the protein
        spectrum:    Spectrum object
    kwargs:
        ion:         str ion type to determine which. Options are 'b', 'y'. Default='b' 
    Ouptut:
        Score object
    '''
    starting_pos = old_entry.starting_position if ion == 'b' else old_entry.starting_position - 1
    ending_pos = old_entry.ending_position + 1 if ion == 'b' else old_entry.ending_position
    if starting_pos < 0 or ending_pos > len(prot) - 1:
        return old_entry

    # check for negative lengths
    if starting_pos > ending_pos or ending_pos < starting_pos:
        return old_entry

    mer_seq = prot[starting_pos:ending_pos+1]
    newmer = Kmer(mer_seq, prot, starting_pos)
    s = Scores(spectrum, newmer)
    s.score()

    return s

def extend_kmer(spectrum: Spectrum, sequence: str, kmer: Scores, ion: str, stall_length=3) -> Scores:
    '''
    Extend a kmer until the score tells us that the adding amino acids doens't make it a better alignment
    
    Inputs:
        spectrum:       Spectrum object
        sequence:       str The full protein sequence we are pulling amino acids from 
        kmer:           Scores object
        ion:            str the ion type we are looking at. Should be 'b' or 'y'
    kwargs:
        stall_length:   int the number of iterations a subsequence is allowed to go witth 
                        no increase in score before finishing kmer growth on a certain kmer. Default=3
    Outputs
        Score object
    '''
    if ion.lower() not in ['b', 'y']:
        return kmer
    score_key = 'b_score' if ion.lower() == 'b' else 'y_score'
    getval = lambda x: getattr(x, score_key) 
    # keep track of the last time a score increased
    last_maintenance = kmer
    # keep going until we run out of extension
    while stall_length > 0:
        updated = new_entry(kmer, sequence, spectrum, ion=ion)
        if getval(updated) > getval(kmer) and getval(updated) > 0 and updated.kmer.n != kmer.kmer.n:
            last_maintenance = updated
        else: 
            stall_length -= 1
        kmer = updated
    return last_maintenance