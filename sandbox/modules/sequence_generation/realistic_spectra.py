import os
import sys
module_path = os.path.abspath(os.path.join('..'))
if module_path not in sys.path:
    sys.path.append(module_path)
module_path = os.path.abspath(os.path.join('../..'))
if module_path not in sys.path:
    sys.path.append(module_path)
    
from src.spectra.gen_spectra import gen_spectrum
from random import randint, uniform
from collections import namedtuple
import numpy as np

RealisticSpectrum = namedtuple('RealisticSpectrum', ['spectrum', 'abundance', 'precursor_mass'])

def gen_realistic_spectra(sequences: list) -> list:
    '''
    Create spectra that look more like real spectra

    Inputs:
        sequences:      (list) sequences of amino acids to generate spectra for
    Ouputs:
        (list of RealisticSpectrum) the spectrum that is closer to a realistic one
    '''
    realistic_spectra = []
    for seq in sequences:
        spec_props = gen_spectrum(seq)
        spec = spec_props['spectrum']
        precursor = spec_props['precursor_mass']
        # Mess with it
        # 1. Drop out peaks
        dropout_rate = randint(60, 85) # rate found from experiments
        dropouts = [randint(0, 100) < dropout_rate for _ in range(len(spec))]
        leftover_peaks = [spec[i] for i in range(len(spec)) if not dropouts[i]]

        # 2. introduce mass errors
        for i in range(len(leftover_peaks)):
            factor = 1 if randint(0, 10) < 5 else -1
            leftover_peaks[i] += factor * np.random.pareto(600) # found from experiments

        # 3. Introduce noise
        leftover_peaks += [uniform(0, max(leftover_peaks) + 100) for _ in range(100-len(leftover_peaks))]

        # 4. pick the abundance
        abundances = np.random.pareto(1, len(leftover_peaks)) * 2000 # found from experiments
        

        realistic_spectra.append(RealisticSpectrum(leftover_peaks, abundances, precursor))
    
    return realistic_spectra