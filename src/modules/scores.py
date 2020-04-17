from src.modules.spectrum import Spectrum

class Scores():
    '''
    Class that holds score information. Contains ion scores for b and y ions

    Inputs:
        n:      int number of top ion scores to hold. Default=3
    '''
    def __init__(self, n=3):
        self.spectrum = None 
        self.n = n
        self.b_scores = None 
        self.y_scores = None

    @property
    def spectrum(self) -> Spectrum:
        return self.spectrum

    @spectrum.setter
    def spectrum(self, spectrum: list, scan_number: int, ms_level: int): 
        self.spectrum = Spectrum(spectrum, scan_number, ms_level)

    