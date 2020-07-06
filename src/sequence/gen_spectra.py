from src.constants import AMINO_ACIDS, WATER_MASS, SINGLY_CHARGED_B_BASE, SINGLY_CHARGED_Y_BASE, DOUBLY_CHARGED_B_BASE, DOUBLY_CHARGED_Y_BASE

def b_ions(sequence: str, charge=None): 
    '''
    Calculate the masses for b_ions
    
    Inputs:
        sequence:  string amino acid sequence to calculate
    kwargs:
        charge:    int charge to calculate. Possible types are {1, 2}. Default is both
    Outputs:
        list of floats
    '''
    masses = []
    length = len(sequence)
    
    if charge is None or charge == 1:
        #b+
        total = SINGLY_CHARGED_B_BASE
        for i in range (0, length):
            total += AMINO_ACIDS[sequence[i]]
            masses.append(total)
            #Since z (the charge) is equal to one, the total here is the m/z
            
    if charge is None or charge == 2:
        #b++
        total = DOUBLY_CHARGED_B_BASE
        for i in range (0, length):
            total += AMINO_ACIDS[sequence[i]]
            masses.append(total/2)
            
    return masses

def y_ions(sequence: str, charge=None): 
    '''
    Calculate the masses for y_ions
    
    Inputs:
        sequence:  string amino acid sequence to calculate
    kwargs:
        charge:    int charge to calculate. Possible types are {1, 2}. Default is both
    Outputs:
        list of floats
    '''
    masses = []
    length = len(sequence)
    
    if charge is None or charge == 1:
        #y+
        total = SINGLY_CHARGED_Y_BASE
        for i in range (0,length):
            total += AMINO_ACIDS[sequence[length-i-1]]
            masses.append(total)
            
    if charge is None or charge == 2:
        #y++
        total = DOUBLY_CHARGED_Y_BASE
        for i in range (0, length):
            total += AMINO_ACIDS[sequence[length-i-1]]
            masses.append(total/2)
            
    return masses

def calc_masses(sequence: str, charge=None, ion=None) -> (list, float):
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
    # sequence = sequence[:-1]
    length = len(sequence)
    total = WATER_MASS
    for i in range(length):
        total +=  AMINO_ACIDS[sequence[i]]

    pre_mz_charge = 2 if charge is None else charge
    pre_mz = (total+pre_mz_charge*1.0072764)/pre_mz_charge   
    
    if ion is None or ion == 'b': 
        masses += b_ions(sequence, charge=charge)
        
    if ion is None or ion == 'y': 
        masses += y_ions(sequence, charge=charge)
        
    return masses, pre_mz

def max_mass(seqeunce: str, ion: str, charge: int) -> float:
    '''
    Calculate the maximum mass of a sequence of an ion type and charge

    Inputs:
        sequence:   (str) the sequence to generate the max mass for
        ion:        (str) the ion type for which we calculate the mass. Options: 'b', 'y'
        charge:     (int) the charge to calculate the mass for. Options: 1, 2
    Outputs:
        (float) the maximum mass
    '''
    # all y ions
    if ion == 'y':
        total = SINGLY_CHARGED_Y_BASE if charge == 1 else DOUBLY_CHARGED_Y_BASE
        total += sum([AMINO_ACIDS[aa] for aa in seqeunce])

        # divide by 2 if doubly charged
        if charge == 2:
            total /= 2
        return total

    # otherwise do the b
    total = SINGLY_CHARGED_B_BASE if charge == 1 else DOUBLY_CHARGED_B_BASE
    total += sum([AMINO_ACIDS[aa] for aa in seqeunce])

    # divide by 2 if doubly
    if charge == 2:
        total /= 2
    
    return total


def gen_spectrum(sequence: str, charge=None, ion=None) -> list:
    '''
    Generate a spectrum for a single sequence. Includes singly and doubly charged masses
    
    Inputs:
        sequence: string amino acid sequence to calculate spectra for
    kwargs:
        charge:   int charge value to calculate masses for. Possible types are {1, 2}. Default is both
        ion:      string ion type to calculate masses for. Possible types are {'b', 'y'}. Default is both
    Outputs:
        dictionary with the following values 
        {
            'spectrum': list of floats,
            'precursor_mass': float,
        }
    '''
    
    this_entry = {}
    masses, pre_mz = calc_masses(sequence, charge=charge, ion=ion)
    this_entry['spectrum'] = masses
    this_entry['precursor_mass'] = pre_mz
    return this_entry

def gen_spectra(sequences: list, charge=None, ion=None) -> list:
    '''
    Generates mass spectra for a list of sequences. Includes singly and doubly charged masses

    Inputs:
        sequences: list of strings sequences to generate spectra for
    kwargs:
        charge:   int charge value to calculate masses for. Possible types are {1, 2}. Default is both
        ion:      string ion type to calculate masses for. Possible types are {'b', 'y'}. Default is both
    Outputs:
        list of dictionaries of the form {'spectrum': list of floats, 'precursor_mass': float}
    '''
    return [gen_spectrum(seq, charge=charge, ion=ion) for seq in sequences]