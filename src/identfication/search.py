from src.types.objects import Spectrum
import math

def ppm_to_da(observed: float, ppm_tolerance: float) -> float:
    '''
    Calculate the absolute boundary value from an observed mass. Value returned is in Da
    
    Inputs:
        observed:      (float) the observed mass used to calculate the tolerances
        ppm_tolerance: (float or int) the tolerance in ppm to use to convert to Da
    Outputs:
        float value in Da 
    '''
    return abs((ppm_tolerance / 1000000)*observed)

def search_kmers_hash(observed: Spectrum, kmers: dict, tolerance: float) -> list:
    '''
    Search through all masses and saved kmers to find masses that are within our tolerance
    
    Inputs:
        spectrum:    (Spectrum) what to sequence
        allbasemers: (dict of list of MassSequence) all of the basemers made from the function 'make_all_base_mers_hash'
        tolerance:   (float) the ppm tolerance to accept for each mass
    Outputs:
        list of MassSequence for all masses that were in the acceptable range of an observed mass
    '''
    hits = []
    for mass in observed.spectrum:
        tol = ppm_to_da(mass, tolerance)
        lb_mass = mass - tol
        ub_mass = mass + tol
        lb_mass_key = math.floor(lb_mass)
        ub_mass_key = math.floor(ub_mass)
        hits += [x.sequence for x in kmers[ub_mass_key] if lb_mass <= x.mass <= ub_mass]
        if lb_mass_key != ub_mass_key:
            hits += [x.sequence for x in kmers[ub_mass_key] if lb_mass <= x.mass <= ub_mass]
            
    return hits
    