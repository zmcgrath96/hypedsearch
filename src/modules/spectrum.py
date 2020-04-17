class Spectrum: 
    '''
    Class that holds spectrum information.
    Should be used as a base class to hold more information

    Inputs:
        spectrum:       list of floats the mass spectrum to hold
        scan_number:    int the scan number in the ms run
        ms_level:       int the ms level to hold
    '''
    def __init__(self, spectrum: list, scan_number: int, ms_level: int):
        self.spectrum = spectrum
        self.scan_number = scan_number
        self.ms_level = ms_level

    @property
    def spectrum(self):
        return self.spectrum

    @property
    def scan_number(self):
        return self.scan_number

    @property 
    def ms_level(self):
        return self.scan_number