from src.modules.spectrum import Spectrum

class Kmer: 
    '''
    Class to hold kmer information and do any transformations

    Public Methods:
        get_spectrum:    get the spectra of this kmer. Can be b or y, singly or doubly, or all of these

    Properties:
        n:              (int) length of sequence
        sequence:       (str) sequence of amino acids
        protein:        (str) name of the protein this sequence is from
        spectrum:       (Spectrum) the Spectrum object for this kmer. Has b and y singly and doubly charged ions. For more fine-grained options
                        call the 'get_spectrum' method with desired parameters
    '''

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

    #This is the mass of water. Adding the mass of water to the sum of all the residue masses gives the mass of the peptide.
    WATER_MASS = 2 * 1.007825035 + 15.99491463 

    SINGLY_CHARGED_Y_BASE = 3 * 1.007825035 + 15.99491463 - 0.0005486 #for the OH to turn the residue CO on the C-terminus into COOH + 1 proton to make NH into NH2 and 1 proton make positively charged
    DOUBLY_CHARGED_Y_BASE = 4 * 1.007825035 + 15.99491463 - 2 * 0.0005486 #another proton to make doubly charged

    SINGLY_CHARGED_B_BASE = 1.007825035 - 0.0005486 #for the H to turn the residue NH on the N-terminus into NH2
    DOUBLY_CHARGED_B_BASE = 2 * 1.007825035 - 2 * 0.0005486 #adding one more proton this time to make it doubly charged

    DEFAULT_ABUNDANCE = 500

    def __init__(self, sequence: str, protein: str, start_position: int):
        self.sequence = sequence
        self.protein = protein
        self.n = len(sequence)
        self.start_position = start_position
        self.end_position = start_position + self.n - 1
        spectrum, precursor_mass = self.__calc_masses()
        self.spectrum = Spectrum(spectrum, 0, 2, precursor_mass, [self.DEFAULT_ABUNDANCE for _ in range(len(spectrum))])

    ############################# Publid Methods #############################
    def get_spectrum(self, charge=None, ion=None) -> list:
        '''
        Generate a spectrum for a single sequence. Includes singly and doubly charged masses
        
        Inputs:
            sequence: string amino acid sequence to calculate spectra for
        kwargs:
            charge:   int charge value to calculate masses for. Possible types are {1, 2}. Default is both
            ion:      string ion type to calculate masses for. Possible types are {'b', 'y'}. Default is both
        Outputs:
            Spectra object
        '''
        
        masses, pre_mz = self.__calc_masses(charge=charge, ion=ion)
        masses.sort()
        return Spectrum(masses, 0, 2, pre_mz, [self.DEFAULT_ABUNDANCE for _ in range(len(masses))])


    ############################# Private Methods #############################
    def __calc_masses(self, charge=None, ion=None) -> (list, float):
        '''
        Calculate the molecular weight (Da) of an Amino Acid sequence
        
        Inputs:
            sequence:   string amino acid sequence to calculate mass of
        kwargs:
            charge:     int charge to calculate. Possible values are {1, 2}. Default is both
            ion:        string ion type to calculate. Possible values are {'b', 'y'}. Default is both
        Output:
            (masses, precursor_mass)
            masses:         list of floats of spectrum calculated
            precursor_mass: float precursor mass of the entire amino acid sequence
        '''
        masses = []

        length = len(self.sequence)
        total = self.WATER_MASS
        for i in range(length):
            total += self.AMINO_ACIDS[self.sequence[i]]

        pre_mz_charge = 2 if charge is None else charge
        pre_mz = (total+pre_mz_charge*1.0072764)/pre_mz_charge   
        
        if ion is None or ion == 'b': 
            masses += self.__b_ions(charge=charge)
            
        if ion is None or ion == 'y': 
            masses += self.__y_ions(charge=charge)
            
        return masses, pre_mz

    # ion generation
    def __b_ions(self, charge=None) -> list: 
        '''
        Calculate the masses for b_ions
        
        kwargs:
            charge:    int charge to calculate. Possible types are {1, 2}. Default is both
        Outputs:
            list of floats
        '''
        masses = []
        length = len(self.sequence)
        
        if charge is None or charge == 1:
            #b+
            total = self.SINGLY_CHARGED_B_BASE
            for i in range (0, length):
                total += self.AMINO_ACIDS[self.sequence[i]]
                masses.append(total)
                #Since z (the charge) is equal to one, the total here is the m/z
                
        if charge is None or charge == 2:
            #b++
            total = self.DOUBLY_CHARGED_B_BASE
            for i in range (0, length):
                total += self.AMINO_ACIDS[self.sequence[i]]
                masses.append(total/2)
                
        return masses

    def __y_ions(self, charge=None) -> list: 
        '''
        Calculate the masses for y_ions
        
        kwargs:
            charge:    int charge to calculate. Possible types are {1, 2}. Default is both
        Outputs:
            list of floats
        '''
        masses = []
        length = len(self.sequence)
        
        if charge is None or charge == 1:
            #y+
            total = self.SINGLY_CHARGED_Y_BASE
            for i in range (0,length):
                total += self.AMINO_ACIDS[self.sequence[length-i-1]]
                masses.append(total)
                
        if charge is None or charge == 2:
            #y++
            total = self.DOUBLY_CHARGED_Y_BASE
            for i in range (0, length):
                total += self.AMINO_ACIDS[self.sequence[length-i-1]]
                masses.append(total/2)
                
        return masses
