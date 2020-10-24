from src.scoring import mass_comparisons
from src.objects import Spectrum
from src.utils import ppm_to_da

from src import gen_spectra

import numpy as np
from bisect import bisect

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

def score_sequence(observed: list, theoretical: list, ppm_tolerance=20) -> float:
    '''
    Score a mass spectrum to a substring of tagged amino acids

    Inputs:
        observed:       (list) the mass spectrum to score
        theoretical:    (list) theoretical spectrum
    kwargs:
        ppm_tolerance:  (int) the tolerance to accepted while scoring. Default=20
    Outputs:
        (float) the score generated from this comparison
    '''

    return mass_comparisons.optimized_compare_masses(observed, theoretical, ppm_tolerance=ppm_tolerance)


def backbone_score(observed: Spectrum, reference: str, ppm_tolerance: int) -> int:
    '''
    Scoring algorithm based on backbone coverage of the reference. The scoring algorithm 
    returns a number between 0 and 100/len(observed) + 3*(len(reference)-1). The calculation of the score is as follows:
    
        1. A percentage is given for the number of bond sites successfully identified
        2. For each bond site that has > 1 ion that describes it, an extra point is awarded. 
        
    Example:
        reference:   ABCDE, 4 junctions to describe
        observed:    ions: b1+, y1++, y2+, b4+
        
        ions for A, E, DE, ABCD found. Coverage = A**DE with DE described by both E and ABCD
        
        Score = %(3/4) + 1 = 75/len(observed) + 1 = 4.75 
        
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
    # extrapoints = sum([jc - 1 if jc > 1 else 0 for jc in jcount])
    
    # make it a function of the number of observed peaks matched
    return jcoverage# + extrapoints

def intensity_backbone_score(observed: Spectrum, reference: str, ppm_tolerance: int) -> int:
    '''
    Scoring algorithm based on backbone coverage of the reference. The scoring algorithm 
    returns a number between 0 and 100/len(observed) + 3*(len(reference)-1). The calculation of the score is as follows:
    
        1. A percentage is given for the number of bond sites successfully identified
        2. For each bond site that has > 1 ion that describes it, an extra point is awarded. 
        
    Example:
        reference:   ABCDE, 4 junctions to describe
        observed:    ions: b1+, y1++, y2+, b4+
        
        ions for A, E, DE, ABCD found. Coverage = A**DE with DE described by both E and ABCD
        
        Score = %(3/4) + 1 = 75/len(observed) + 1 = 4.75 
        
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

    # keep track of the abundances that contribute to our score
    ided_abundances = 0
    
    for ion in ['b', 'y']:
        for charge in [1, 2]:
            # take off the trailing or leading amino acid from the reference according to ion type
            singled_seq = reference[:-1] if ion == 'b' else reference[1:]

            # get the m/z peaks
            peaks = gen_spectra.gen_spectrum(singled_seq, charge=charge, ion=ion)['spectrum']
            
            # go through each peak and try and match it to an observed one
            for i in range(len(peaks)):

                # take tolerance into account 
                da_tol = ppm_to_da(peaks[i], ppm_tolerance)

                # get hits
                peak_hits = list(map(lambda idx_x: idx_x[0] if peaks[i] - da_tol <= idx_x[1] <= peaks[i] + da_tol else None, enumerate(observed.spectrum)))

                # remove None from peak hits
                peak_hits = [idx for idx in peak_hits if idx is not None]

                # if the len > 1, increment
                if len(peak_hits) > 0:
                    jcount[i] += 1
                    ided_abundances += sum([observed.abundance[idx] for idx in peak_hits])
    
    jcoverage = sum([1 if jc > 0 else 0 for jc in jcount])
    #extrapoints = sum([jc - 1 if jc > 1 else 0 for jc in jcount]) // 3
    
    # make it a function of the number of observed peaks matched
    #jcoverage /= len(observed.spectrum)
    return jcoverage #+ extrapoints) * (ided_abundances / observed.total_intensity)


def ion_backbone_score(observed: Spectrum, reference: str, ion: str, ppm_tolerance: int) -> float:
    '''
    Scoring algorithm based on backbone coverage of the reference. The scoring algorithm 
    returns a number between 0 and (100/+ 3*(len(reference)-1))/len(observed) . The calculation of the score is as follows:
    
        1. A percentage is given for the number of bond sites successfully identified
        2. For each bond site that has > 1 ion that describes it, an extra point is awarded. 
        
    Example:
        reference:   ABCDE, 4 junctions to describe
        observed:    ions: b1+, y1++, y2+, b4+
        
        ions for A, E, DE, ABCD found. Coverage = A**DE with DE described by both E and ABCD
        
        Score = %(3/4) + 1 = 75/len(observed) + 1 = 4.75 
        
    Inputs:
        observed:       (Spectrum) spectrum being scored against
        reference:      (str) reference amino acid sequence being scored against the spectrum
        ion:            (str) the ion type to focus on. Options are 'b' or 'y'
        ppm_tolerance:  (int) tolerance to allow in ppm for each peak
    Outputs:
        (float) score according to the formula 
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
    divider = 1 if len(jcount) <= 1 else len(jcount)
    jcoverage = int(100 * sum([1 if jc > 0 else 0 for jc in jcount]) / divider)
    #extrapoints = sum([jc - 1 if jc > 1 else 0 for jc in jcount])
    
    return jcoverage #+ extrapoints

def intensity_ion_backbone_score(observed: Spectrum, reference: str, ion: str, ppm_tolerance: int) -> float:
    '''
    Scoring algorithm that factors in how much of the backbone is identified along with the abundance each peak 
    contributes to the total score. The formula is:

        count(# backbone cleavages found) + 1 for every additional hit of an already identifed cleavage * total percentage
        of intensity covered by identfied peaks

    Example:
        reference:   ABCDE, 4 junctions to describe
        observed:    ions: b1+, y1++, y2+, b4+, with relative intensities of (100, 200, 100, 50) of a total 600
        
        ions for A, E, DE, ABCD found. Coverage = A**DE with DE described by both E and ABCD
        
        Score = %(3/4) + 1 = 75 + 1 = 76 *(100 + 200 + 100 + 50)/600 = 57 

    Inputs:
        observed:       (Spectrum) spectrum being scored against
        reference:      (str) reference amino acid sequence being scored against the spectrum
        ion:            (str) the ion type to focus on. Options are 'b' or 'y'
        ppm_tolerance:  (int) tolerance to allow in ppm for each peak
    Outputs:
        (float) score according to the formula 
    '''
    # check to see if observed is nothing
    if len(observed.spectrum) == 0:
        return 0

    # keep track of the junction (bond) sites found
    jcount = [0 for _ in range(len(reference)-1)]

    # keep track of the abundances that contribute to our score
    ided_abundances = 0
    
    for charge in [1, 2]:

        # take off the trailing or leading amino acid from the reference according to ion type
        singled_seq = reference[:-1] if ion == 'b' else reference[1:]

        # get the m/z peaks
        peaks = gen_spectra.gen_spectrum(singled_seq, charge=charge, ion=ion)['spectrum']
        
        # go through each peak and try and match it to an observed one
        for i in range(len(peaks)):

            # take tolerance into account 
            da_tol = ppm_to_da(peaks[i], ppm_tolerance)

            # get hits
            peak_hits = list(map(lambda idx_x: idx_x[0] if peaks[i] - da_tol <= idx_x[1] <= peaks[i] + da_tol else None, enumerate(observed.spectrum)))

            # remove None from peak hits
            peak_hits = [idx for idx in peak_hits if idx is not None]

            # if the len > 1, increment
            if len(peak_hits) > 0:
                jcount[i] += 1
                ided_abundances += sum([observed.abundance[idx] for idx in peak_hits])

    # if any entry has at least 1, that bond has been identified at least once
    jcoverage = sum([1 if jc > 0 else 0 for jc in jcount])

    # if an entry has more than 1, we give it extra points
    #extrapoints = sum([jc - 1 if jc > 1 else 0 for jc in jcount]) // 3
    
    return jcoverage# + extrapoints) * (ided_abundances / observed.total_intensity)

def ion_intensity_percentage(observed: Spectrum, reference: str, ppm_tolerance: int, ion: str) -> int:
    '''
    The percentage the matched m/z values make up of the overall intensity. For example,
    if 4 peaks are matched between the observed and reference and those peaks make up 10% of the
    total intensity, 10 is returned. Spectrum m/z values should be sorted low to high with the intensity
    values indexed in the same order as the m/z values

    Inputs:
        observed:       (Spectrum) observed spectrum
        reference:      (str) the reference sequence
        ppm_tolerance:  (int) the tolerance (in ppm) to allow when matching peaks
    Output:
        (int) rounded percentage value
    '''
    if len(reference) == 0:
        return 0

    # generate the spectrum for the reference sequence
    refspec = gen_spectra.gen_spectrum(reference, ion=ion)['spectrum']
        
    def boundaries(mass):
        tol = ppm_to_da(mass, ppm_tolerance)
        return [mass - tol, mass + tol]
                
    # calculate the boundaries for each of the reference masses for binary search
    reference_boundaries = []
    for mz in refspec:
        reference_boundaries += boundaries(mz)
        
    # get all indices where the observed gets a hit
    idxes = [i for i, mz in enumerate(observed.spectrum) if bisect(reference_boundaries, mz) % 2]

    # take the intensities at idxes and sum them
    score = int( 100 * sum([observed.abundance[i] for i in idxes]) / observed.total_intensity)
    
    return score


def hybrid_score(observed: Spectrum, hybrid_seq: str, ppm_tolerance: int) -> float:
    '''
    A stricter score for hybrid sequences. Ions found on the relative side of the hybrid junction
    (left for b, right for y) will only be awarded .5 points, where as an ion found on the other 
    side of the junction (right for b, left for y) will be awarded a full point

    Example:
        ABC-DEF

        b-ions: A, C, E
        y-ions: D, A

        Score .5(bA) + .5(bC) + 1(bE) + .5 (yD) + 1(yA) = 3.5

    Inputs:
        observed:       (Spectrum) observed spectrum
        hybrid_seq:     (str) the hybrid sequence to score against
        ppm_tolerance:  (int) the tolerance to allow when searching for a mass
    Outputs:
        (float) the score
    '''
    if '-' not in hybrid_seq and '(' not in hybrid_seq and ')' not in hybrid_seq:
        return 0

    score = 0

    # get a non hybrid for scoring purposes
    non_hyb = hybrid_seq.replace('-', '').replace('(', '').replace(')', '')

    # get the index to the left of which b ions will only get .5 points
    b_split = hybrid_seq.index('-') if '-' in hybrid_seq else hybrid_seq.index('(')

    # get the index to the right of which y ions will only get .5 points
    y_split = len(hybrid_seq) - hybrid_seq.index('-') if '-' in hybrid_seq else len(hybrid_seq) - hybrid_seq.index(')') 

    # generate b and y separately to be sure
    b_spec = sorted(gen_spectra.gen_spectrum(non_hyb, ion='b')['spectrum'])
    y_spec = sorted(gen_spectra.gen_spectrum(non_hyb, ion='y')['spectrum'])

    # convert the spectra into lists of tuples
    gen_range = lambda x: (x - ppm_to_da(x, ppm_tolerance), x + ppm_to_da(x, ppm_tolerance))
    b_ranges = [gen_range(x) for x in b_spec]
    y_ranges = [gen_range(x) for x in y_spec]

    # do a merge search where we linearly search each mass in the observed twice
    b_range_i = 0
    observed_i = 0
    while b_range_i < len(b_ranges) and observed_i < len(observed.spectrum):

        # if observed is larger than the range, increment range
        if observed.spectrum[observed_i] > b_ranges[b_range_i][1]:
            b_range_i += 1

        # if observed is smaller than the range, increment observed
        elif observed.spectrum[observed_i] < b_ranges[b_range_i][0]:
            observed_i += 1

        # otherwise its in the range, see what to increment score by, and increment observed
        else:
            score += 1 if b_range_i >= b_split else .5
            observed_i += 1

    y_range_i = 0
    observed_i = 0
    while y_range_i < len(y_ranges) and observed_i < len(observed.spectrum):

        # if observed is larger than the range, increment range
        if observed.spectrum[observed_i] > y_ranges[y_range_i][1]:
            y_range_i += 1

        # if observed is smaller than the range, increment observed
        elif observed.spectrum[observed_i] < y_ranges[y_range_i][0]:
            observed_i += 1

        # otherwise its in the range, see what to increment score by, and increment observed
        else:
            score += 1 if y_range_i >= y_split else .5
            observed_i += 1

    return score

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