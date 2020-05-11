class Spectrum: 

    AMINO_ACIDS = {
        "A": 71.037114,
        "R": 156.101111,
        "N": 114.042927,
        "D": 115.026943,
        "C": 103.009185,
        "E": 129.042593,
        "Q": 128.058578,
        "G": 57.021464,
        "H": 137.058912,
        "I": 113.084064,
        "L": 113.084064,
        "K": 128.094963,
        "M": 131.040485,
        "F": 147.068414,
        "P": 97.052764,
        "S": 87.032028,
        "T": 101.047679,
        "U": 150.95363,
        "W": 186.079313,
        "Y": 163.06332,
        "V": 99.068414,
        "X": 0, # added to ignore. TODO: figure out what to do with it
        "B": 0, # added to ignore. TODO: figure out what to do with it
        "Z": 0, # added to ignore. TODO: figure out what to do with it
    }

    '''
    Class that holds spectrum information.
    Should be used as a base class to hold more information

    Inputs:
        spectrum:       list of floats the mass spectrum to hold
        scan_number:    int the scan number in the ms run
        ms_level:       int the ms level to hold
    kwargs:
        abundance:      list of floats of the corresponding abundance of each peak
    '''
    def __init__(self, spectrum: list, scan_number: int, ms_level: int, abundance=[]):
        self.spectrum = spectrum
        self.abundance = abundance
        self.scan_number = scan_number
        self.ms_level = ms_level

    ################################## PUBLIC METHODS ##################################

    def score(self, sequence: str):
        '''
        Score my spectrum to a substring of tagged amino acids

        Inputs:
            pepspec:    list of floats the mass spectrum to score
            subseq:     str of tagged amino acids to score the spectrum against
        Outputs:
            (b_score, y_score): (float, float) the b and y ion scores generated from this comparison
        '''
        kmerspec_b = gen_spectra.gen_spectrum(subseq, ion='b')['spectrum']
        kmerspec_y = gen_spectra.gen_spectrum(subseq, ion='y')['spectrum']
        b_score = self.__compare_masses(self.spectrum, kmerspec_b)
        y_score = self.__compare_masses(self.spectrum, kmerspec_y)
        return (b_score, y_score)

    ################################## PRIVATE METHODS ##################################

    def __compare_masses(spectrum: list, reference: list) -> float:
        '''
        CREATED APRIL 27 2020
        Score two spectra against eachother. Simple additive scoring with bonuses for streaks
        Divides by the length of the reference to make it length biased for the reference

        Note:   the difference between this one and the other April one is this one divides 
                by the length of the spectrum. This is because all extended kmers will 
                be getting longer, so we want the score to increase with the length, not stay the same

        Inputs:
            spectrum:   list of floats (from mass spectra)
            reference:  list of floats (calculated from protein sequence)
        Outputs:
            score:      float score 
        '''
        if len(spectrum) == 0 or len(reference) == 0:
            return 0.0
        streak = 0
        last = True
        score = 0
        max_streak = 0
        for refmass in reference:
        
            if refmass in spectrum:
                if last == True:
                    streak += 1
                    max_streak = max([streak, max_streak])
                score += 1
                last = True 

            else:
                streak = 0
                last = False
        
        score += max_streak
        score /= float(len(spectrum))
        return score 