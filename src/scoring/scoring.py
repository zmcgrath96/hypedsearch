from src.scoring import mass_comparisons
from src.spectra import gen_spectra
from src.types.objects import Spectrum
from src.utils import ppm_to_da

import numpy as np

def score_subsequence(pepspec: list, subseq: str, ppm_tolerance=20) -> (float, float):
    '''
    Score a mass spectrum to a substring of tagged amino acids

    Inputs:
        pepspec:        (list of floats) the mass spectrum to score
        subseq:         (str) amino acids to score the spectrum against
    kwargs:
        ppm_tolerance:  (int) the tolerance to accepted while scoring. Default=20
    Outputs:
        (b_score, y_score): (float, float) the b and y ion scores generated from this comparison
    '''
    kmerspec_b = gen_spectra.gen_spectrum(subseq, ion='b')['spectrum']
    kmerspec_y = gen_spectra.gen_spectrum(subseq, ion='y')['spectrum']
    b_score = mass_comparisons.optimized_compare_masses(pepspec, kmerspec_b, ppm_tolerance=ppm_tolerance)
    y_score = mass_comparisons.optimized_compare_masses(pepspec, kmerspec_y, ppm_tolerance=ppm_tolerance)
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
    if len(reference) < 2:
        return 0

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
    
    # make it a function of the number of observed peaks matched
    jcoverage /= len(observed.spectrum)
    return jcoverage + extrapoints


def ion_backbone_score(observed: Spectrum, reference: str, ion: str, ppm_tolerance: int) -> float:
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
        ion:            (str) the ion type to focus on. Options are 'b' or 'y'
        ppm_tolerance:  (int) tolerance to allow in ppm for each peak
    '''
    jcount = [0 for _ in range(len(reference)-1)]
    
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
    
    # make it a fuction of the number of observed peaks matched
    jcoverage /= len(observed.spectrum)
    return jcoverage + extrapoints

def precursor_distance(observed_precursor: float, reference_precursor: float) -> float:
    '''
    The absolute distance between two precursor masses

    Inputs:
        observed_precursor:     (float) observed precursor mass
        reference_precursor:    (float) reference precursor mass
    Outputs:
        (float) absolute distance between the two
    '''
    return abs(observed_precursor - reference_precursor)


def xcorr(observed: np.ndarray, reference: np.ndarray, value=50) -> float:
    '''
    An xcorr scoring algorithm 
    
    Formula: dot(reference, y'), y' = observed - (sum(-75, 75)observed[i]/150)
    
    Inputs should be sparsely populated np arrays. Indexing should observe the formula
    
    idx = int(m/w), m is mass, w is bin width (aka tolerance)
    
    Inputs:
        observed:    (np.ndarray) list of peaks normalized to value
        reference:   (np.ndarray) list of peaks normalized to value
    kwargs: 
        value:       (number) value given to the normalized peaks. Default=50
    Outputs:
        (float) the score of the formula shown above
    '''
    def sum_term(i):
        min_idx = max(0, i-75)
        max_idx = min(len(observed) - 1, i + 75)
        return np.sum(observed[min_idx:max_idx])/(max_idx - min_idx)

    # calculate the y prime term 
    y_prime = np.asarray([observed[i] - sum_term(i) for i in range(len(observed))])
        
    # fill y prime or reference to be the same size. fill with zeros
    if len(y_prime) < len(reference):
        y_prime = np.concatenate((y_prime, [0 for _ in range(len(reference) - len(y_prime))]))
    else:
        reference = np.concatenate((reference, [0 for _ in range(len(y_prime) - len(reference))]))
        
    return np.dot(reference, y_prime)/(sum(reference)*value)


def overlap_score(l: list, ion: str) -> float:
    '''
    A score based on how well the sequences that make up the list overlap

    Example:   
        l: [ABC, ABCDE, ABCDLMNO, ABCDEF]

        + 3 for ABC overlapping all of them
        + 1 for ABCDE overlapping ABCDEF
        - .5 for ABCDEF not overlapping ABCDLMNO
        Final score: 3.5

    Inputs:
        l:      (list) sequences checking for overlap
        ion:    (str) ion type. If y, strings compared right to left
    Outputs:
        (float) score of the overlap
    '''
    l.sort(key=lambda x: len(x))
    
    score = 0
    
    for i, e in enumerate(l[:-1]):
        
        # if its y, we want to go from right to left
        if ion == 'y':
            e = e[::-1]
        
        for j, e2 in enumerate(l[i+1:]):
            
            # if its y, we want to go from right to left
            if ion == 'y':
                e2 = e2[::-1]
            
            score += 1 if e2[:len(e)] == e else -.5
            
    return score