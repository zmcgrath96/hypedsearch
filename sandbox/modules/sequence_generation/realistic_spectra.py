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

def print_missing_peaks(leftovers, sequence):
    bs = gen_spectrum(sequence, charge=1, ion='b')['spectrum']
    bd = gen_spectrum(sequence, charge=2, ion='b')['spectrum']
    ys = gen_spectrum(sequence, charge=1, ion='y')['spectrum']
    yd = gen_spectrum(sequence, charge=2, ion='y')['spectrum']

    bsp = [i for i in range(len(bs)) if bs[i] in leftovers]
    bdp = [i for i in range(len(bd)) if bd[i] in leftovers]
    ysp = [i for i in range(len(ys)) if ys[i] in leftovers]
    ydp = [i for i in range(len(yd)) if yd[i] in leftovers]

    rbs = len([x for x in bs if x in leftovers])
    rbd = len([x for x in bd if x in leftovers])
    rys = len([x for x in ys if x in leftovers])
    ryd = len([x for x in yd if x in leftovers])

    print(f'Count of remaining ions by type:\nb+: {rbs} \t b++: {rbd} \t y+: {rys} \t y++: {ryd}')
    print(f'Remaing ion types by position:\nb+: {bsp}\nb++: {bdp}\ny+: {ysp}\ny++: {ydp}')

def gen_realistic_spectra(sequences: list, DEBUG=False) -> list:
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
        DEBUG and print(f'Dropout rate: {dropout_rate}')
        dropouts = [randint(0, 100) < dropout_rate for _ in range(len(spec))]
        leftover_peaks = [spec[i] for i in range(len(spec)) if not dropouts[i]]
        DEBUG and print(f'Number of peaks remaining: {len(leftover_peaks)}/{len(spec)}')
        DEBUG and print_missing_peaks(leftover_peaks, seq)
        # 2. introduce mass errors
        for i in range(len(leftover_peaks)):
            factor = 1 if randint(0, 10) < 5 else -1
            error = factor * np.random.pareto(600)
            leftover_peaks[i] += error  #found from experiments

        # 3. Introduce noise
        leftover_peaks += [uniform(0, max(leftover_peaks) + 100) for _ in range(100-len(leftover_peaks))]

        # 4. pick the abundance
        abundances = list(np.random.pareto(1, len(leftover_peaks)) * 2000) # found from experiments
        
        realistic_spectra.append(RealisticSpectrum(leftover_peaks, abundances, precursor))
    
    return realistic_spectra