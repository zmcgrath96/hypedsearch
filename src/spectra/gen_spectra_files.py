import argparse 
import json
from spectra import gen_spectra, write_spectra
from sequences.gen_k_mers import k_mers
from utils.utils import make_valid_dir_string, make_dir

'''generate 

DESC:
    generate kmers and the spectra files associated with them 
Inputs:
    sequences: list of dictionaries of the form {'name': str, 'sequence': str}. Sequences to generate
                kmers from
    window_sizes: list of ints size of kmers to generate
kwargs:
    save_dir: str the directory in which to save all the spectra files. Default=./
    compress: bool compress the spectra files. Default=True
Outputs:
    list of strs of the file names/paths generated
'''
def generate(sequences, window_sizes, save_dir='./', compress=True):
    output_files = []
    save_dir = make_valid_dir_string(save_dir) + 'spectra/'
    make_dir(save_dir)
    
    for window_size in window_sizes:
        print('Generating {}-mer spectra for all proteins...'.format(window_size))
        seq_c = 0
        n_s = len(sequences)
        for sequence in sequences:
            print('generating spectra for sequence {}/{}[{}%]\r'.format(seq_c, n_s, int((float(seq_c)/(float(n_s))) * 100)), end="")
            name = '{}_{}'.format(sequence['name'], window_size)
            kmers = k_mers(sequence['sequence'], window_size)
            spectra = gen_spectra.gen_spectra(kmers)
            output_files.append(write_spectra.write_mzml(name, spectra, output_dir=save_dir, compress=compress))
            seq_c += 1

    return list(set(output_files))
