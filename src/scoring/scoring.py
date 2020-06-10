from src.scoring import mass_comparisons
from src.spectra import gen_spectra
from src.types.objects import Spectrum
from src.utils import ppm_to_da

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


def backbone_score(observed: Spectrum, reference: str, ppm_tolerance: int) -> int:
    '''
    Scoring algorithm based on backbone coverage of the reference. The scoring algorithm 
    returns a number between 0 and 100 + 3*(len(reference)-1). The calculation of the score is as follows:
    
        1. A percentage is given for the number of bond sites successfully identified
        2. For each bond site that has > 1 ion that describes it, an extra point is awarded. 
        
    Example:
        reference:   ABCDE, 4 junctions to describe
        observed:    ions: b1+, y1++, y2+, b4+
        
        ions for A, E, DE, ABCD found. Coverage = A**DE with DE described by both E and ABCD
        
        Score = %(3/4) + 1 = 75 + 1 = 76 
        
    Inputs:
        observed:       (Spectrum) spectrum being scored against
        reference:      (str) reference amino acid sequence being scored against the spectrum
        ppm_tolerance:  (int) tolerance to allow in ppm for each peak
    Outputs:
        (int) score according the the function described above
    '''
    jcount = [0 for _ in range(len(reference)-1)]
    
    for ion in ['b', 'y']:
        for charge in [1, 2]:
            singled_seq = reference[:-1] if ion == 'b' else reference[1:]
            peaks = gen_spectra.gen_spectrum(singled_seq, charge=charge, ion=ion)['spectrum']
            peaks = peaks if ion == 'b' else peaks[::-1]
            for i in range(len(peaks)):
                da_tol = ppm_to_da(peaks[i], ppm_tolerance)
                if any([peaks[i] - da_tol <= obs_peak <= peaks[i] + da_tol for obs_peak in observed.spectrum]):
                    jcount[i] += 1
    
    jcoverage = int(100 * sum([1 if jc > 0 else 0 for jc in jcount]) / len(jcount))
    extrapoints = sum([jc - 1 if jc > 1 else 0 for jc in jcount])
    
    return jcoverage + extrapoints