from src.scoring import mass_comparisons
from src.spectra import gen_spectra

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
    b_score = mass_comparisons.optimized_compare_masses(pepspec, kmerspec_b)
    y_score = mass_comparisons.optimized_compare_masses(pepspec, kmerspec_y)
    return (b_score, y_score)
