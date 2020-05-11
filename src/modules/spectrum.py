class Spectrum: 

    '''
    Class that holds spectrum information.
    Should be used as a base class to hold more information

    Public Methods:


    Properties:
        spectrum:       (list) of floats the mass spectrum to hold
        scan_number:    (int) the scan number in the ms run
        ms_level:       (int) the ms level to hold
        precursor_mass: (float) the precursor mass of the spectrum
        abundance:      (list) of floats of the corresponding abundance of each peak
    '''
    def __init__(self, spectrum: list, scan_number: int, ms_level: int, precursor_mass: int, abundance=[]):
        self.spectrum = spectrum
        self.abundance = abundance
        self.scan_number = scan_number
        self.ms_level = ms_level
        self.precursor_mass = precursor_mass

    ################################## Public Methods ##################################

    ################################## Private Methods ##################################
