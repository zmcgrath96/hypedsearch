from src.constants import AMINO_ACIDS, WATER_MASS, SINGLY_CHARGED_B_BASE, SINGLY_CHARGED_Y_BASE, DOUBLY_CHARGED_B_BASE, DOUBLY_CHARGED_Y_BASE, INTEGER_ORDERED_AMINO_ACIDS

import numpy as np

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


def get_precursor(sequence: str, charge=2) -> float:
    '''
    Calculate JUST the precursor mass of the input sequence at the charge provided.

    Inputs:
        sequence:   (str) the AA sequence to calculate the precursor for
    kwargs:
        charge:     (int) the charge for which to calculate the precursor. Default=2
    Outputs:
        (float) the precursor mass of the sequence
    '''
    total = WATER_MASS
    for aa in sequence:
        total +=  AMINO_ACIDS[aa]

    pre_mz_charge = 2 if charge is None else charge
    return (total+pre_mz_charge*1.0072764)/pre_mz_charge  


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

def gen_min_ordering(sequence: str) -> list:
    '''
    Generates an np array the length of the sequence that is the minimal representation
    of a spectrum (for ordering purposes). Each amino acid is represented as an integer
    (an 8 bit integer by NumPy). The integer is the sorted value (lowest to highest) by mass. 
    For example, G has the lowest mass, so its index is 0. W, the heaviest, has an index of 19.
    This is done for the smallest (memory) representation for sorting a list of k-mers
    (needed for simple DAWG construction)

    Inputs:
        sequence:   (str) the sequence of amino acids to converte
    Outputs:
        (list) a list of 16 bit integers
    '''

    middle = [np.int8(INTEGER_ORDERED_AMINO_ACIDS[aa]) for aa in sequence if aa in INTEGER_ORDERED_AMINO_ACIDS]
    if len(middle) == 0:
        return []
    return middle

   