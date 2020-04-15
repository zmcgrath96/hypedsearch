from src.alignment import mass_comparisons
from src.spectra import gen_spectra

def confidence(b_entry: dict, y_entry: dict) -> float:
    '''
    CREATED 15 APRIL 2020
    Calculate a confidence score from some alignment. The calculation works as follows:
        overlapped_b_ion:   and identified b_ion that are in the identified y_ion range
        overlapped_y_ion:   and identified y_ion that are in the identified b_ion range

        ((# overlapped_b_ions * b_score) + (# overlapped_y_ions * y_score) + # total ions identified) - (# non overlapped_b_ions + # non overlapped_y_ions)

    Inputs:
        (b_entry, y_entry):     dict of the form 
            {starting_position: int, ending_position: int, k: int, b_score: float, y_score: float`}
    Outputs:
        float the score calculated by the above formula
    '''
    b_range = range(b_entry['starting_position'], b_entry['ending_position'] + 1)
    y_range = range(y_entry['starting_position'], y_entry['ending_position'] + 1)
    overlapped_b_ions = sum([1 if x in y_range else 0 for x in b_range])
    overlapped_y_ions = sum([1 if x in b_range else 0 for x in y_range])
    non_overlapped_b_ions = sum([1 if x not in y_range else 0 for x in b_range])
    non_overlapped_y_ions = sum([1 if x not in b_range else 0 for x in y_range])
    
    return (overlapped_b_ions * b_entry['b_score'] + overlapped_y_ions * y_entry['y_score'] + len(b_range) + len(y_range)) - (non_overlapped_b_ions + non_overlapped_y_ions)

def confidence_simple(b_entry: dict, y_entry: dict) -> float:
    '''
    CREATED 15 APRIL 2020
    Calculate a simple confidence score. The calculation is a follows:
        
        score = # overlaps (number of amino acids described by both the b and y ions)
        perfect_score = max(len(b_sequence), len(y_sequence))
        confidence = (score / perfect_score) * 100 

    Inputs:
        (b_entry, y_entry):     dict of the form 
            {starting_position: int, ending_position: int, k: int, b_score: float, y_score: float`}
    Outputs:
        float the score calculated by the above formula
    '''
    b_range = range(b_entry['starting_position'], b_entry['ending_position'] + 1)
    y_range = range(y_entry['starting_position'], y_entry['ending_position'] + 1)
    perfect_score = max(len(b_range), len(y_range))
    score = sum([1 if x in y_range else 0 for x in b_range])
    return (float(score) / float(perfect_score)) * 100

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
    b_score = mass_comparisons.compare_masses(pepspec, kmerspec_b)
    y_score = mass_comparisons.compare_masses(pepspec, kmerspec_y)
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