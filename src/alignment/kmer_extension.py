from src.scoring import mass_comparisons
from src.spectra import gen_spectra
from src.types.objects import ScoredKmer, Spectrum, Kmer
from src.types.database import Entry

def new_entry(old_entry: ScoredKmer, protein: Entry, spectrum: Spectrum, ion='b') -> ScoredKmer:
    '''
    Generate a new entry from the old entry
    
    Input:
        old_entry:   dict entry with k, sequence, b and y scores, start and end positions
        prot:        Entry class instance 
        spectrum:    list spectrum to score against
    kwargs:
        ion:         str ion type to determine which. Options are 'b', 'y'. Default='b' 
    Ouptut:
        Score object
    '''
    starting_pos = old_entry.kmer.start_position if ion == 'b' else old_entry.kmer.start_position - 1
    ending_pos = old_entry.kmer.end_position + 1 if ion == 'b' else old_entry.kmer.end_position
    if starting_pos < 0 or ending_pos > len(protein.sequence) - 1:
        return old_entry

    # check for negative lengths
    if starting_pos > ending_pos or ending_pos < starting_pos:
        return old_entry

    mer_seq = protein.sequence[starting_pos:ending_pos+1]
    mer_spec_b = gen_spectra.gen_spectrum(mer_seq, ion='b')['spectrum']
    mer_spec_y = gen_spectra.gen_spectrum(mer_seq, ion='y')['spectrum']
    longer_kmer = Kmer(old_entry.kmer.k + 1, mer_seq, protein.name, starting_pos, ending_pos)
    new_sk = ScoredKmer( mass_comparisons.compare_masses(spectrum.spectrum, mer_spec_b), mass_comparisons.compare_masses(spectrum.spectrum, mer_spec_y), longer_kmer)
    return new_sk

def extend_kmer(spectrum: Spectrum, protein: Entry, kmer: ScoredKmer, ion: str, stall_length=3) -> dict:
    '''
    Extend a kmer until the score tells us that the adding amino acids doens't make it a better alignment
    
    Inputs:
        spectrum:           Spectrum namedtuple instance
        protein_sequence:   Entry class instance
        kmer:               ScoredKmer namedtuple instance
        ion:                str the ion type we are looking at. Should be 'b' or 'y'
    kwargs:
        stall_length:   int the number of iterations a subsequence is allowed to go witth 
                        no increase in score before finishing kmer growth on a certain kmer. Default=3
    Outputs
        Score object
    '''
    if ion.lower() not in ['b', 'y']:
        return kmer
    score_key = 'b_score' if ion.lower() == 'b' else 'y_score'
    # keep track of the last time a score increased
    last_maintenance = kmer
    # keep going until we run out of extension
    while stall_length > 0:
        updated = new_entry(kmer, protein, spectrum, ion=ion)
        if getattr(updated, score_key) > getattr(kmer, score_key) and getattr(updated, score_key) > 0 and updated.kmer.k != kmer.kmer.k:
            last_maintenance = updated
        else: 
            stall_length -= 1
        kmer = updated
    return last_maintenance