from src.modules.spectrum import Spectrum
from src.modules.kmer import Kmer
class Scores():
    '''
    Class that holds score information. Contains ion scores for b and y ions

    Public Methods:
        score:          function called in order to score the kmer against the spectrum

    Properties: 
        spectrum:       (Spectrum) Spectrum object used in scoring
        kmer:           (Kmer) Kmer object used to score against spectrum
        b_score:        (float) the b ion score. Result of the score method
        y_score:        (float) the y ion score. Result of the score method

    '''
    def __init__(self, specturm: Spectrum, kmer: Kmer):
        self.kmer = kmer
        self.spectrum = specturm
        self.b_score = None 
        self.y_score = None

    ################################## Public Methods ##################################
    def score(self) -> None:
        '''
        Score the kmer property against the spectrum
        '''
        kmerspec_b = self.kmer.get_spectrum(ion='b')
        kmerspec_y = self.kmer.get_spectrum(ion='y')
        self.b_score = self.__compare_masses(kmerspec_b.spectrum)
        self.y_score = self.__compare_masses(kmerspec_y.spectrum)

    ################################## Private Methods ##################################
    def __compare_masses(self, reference: list) -> float:
        '''
        CREATED APRIL 27 2020
        Score two spectra against eachother. Simple additive scoring with bonuses for streaks
        Divides by the length of the reference to make it length biased for the reference

        Note:   the difference between this one and the other April one is this one divides 
                by the length of the spectrum. This is because all extended kmers will 
                be getting longer, so we want the score to increase with the length, not stay the same

        Inputs:
            reference:  list of floats (calculated from protein sequence)
        Outputs:
            score:      float score 
        '''
        if len(reference) == 0:
            return 0.0
        streak = 0
        last = True
        score = 0
        max_streak = 0
        for refmass in reference:
        
            if refmass in self.spectrum:
                if last == True:
                    streak += 1
                    max_streak = max([streak, max_streak])
                score += 1
                last = True 

            else:
                streak = 0
                last = False
        
        score += max_streak
        score /= float(len(self.spectrum))
        return score 