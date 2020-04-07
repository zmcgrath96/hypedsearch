'''
id_spectra.py

Author: Zachary McGrath
Date: 6 April 2020

Find protein and subsequence for every spectra passed in
'''
from file_io import fasta, mzML

def id_spectrum(spectrum: dict, database: list) -> dict:
    '''
    Run a scoring and alignment on a specific spectrum 

    Inputs:
        spectrum:   dictionary with the following attributes
            scan_no:    int scan number
            spectrum:   list of floats
            level:      int ms level
        database:   list of proteins loaded from a .fasta file
    Outputs:
        dict of the following form
        {
            'spectrum': list,
            'alignments': list of dicts,
        }
    '''
    pass

def id_spectra(spectra_files: list, database_file: str) -> dict:
    '''
    Run a scoring and alignment on each spectra passed in

    Inputs:
        spectra_files:  list of strings of file names of spectra
        database_file:  string full path to a .fasta database
    Outputs:
        dict containing the results. All information is keyed by the spectrum file name with scan number appended
    '''
    database = fasta.read(database_file, True)
    results = {}
    for i, spectrum_file in enumerate(spectra_files):
        print('Spectra file {}/{}[{}%]'.format(i, len(spectra_files), int(float(i)/float(len(spectra_files)) * 100)))
        spectra = mzML.read(spectrum_file)
        for spectrum in spectra:
            entry = id_spectra(spectrum, database)
            entry_name = '{}_{}'.format(spectrum_file, spectrum['scan_no'])
            results[entry_name] = entry