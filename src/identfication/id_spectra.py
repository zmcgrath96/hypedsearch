'''
id_spectra.py

Author: Zachary McGrath
Date: 6 April 2020

Find protein and subsequence for every spectra passed in
'''
from src.file_io import fasta, mzML
from src.alignment import search, aligners
from math import ceil

############## Constants ##############
TOP_N = 3
#######################################

def id_spectrum(spectrum: dict, database: list, n=TOP_N) -> dict:
    '''
    Run a scoring and alignment on a specific spectrum 

    Inputs:
        spectrum:   dictionary with the following attributes
            scan_no:    int scan number
            spectrum:   list of floats
            level:      int ms level
        database:   list of proteins loaded from a .fasta file
    kwargs:
        n:          number of alignments to return for each spectrum. Default=3
    Outputs:
        dict of the following form
        {
            'spectrum': list,
            'alignments': list of dicts,
            'scan_no': int
        }
    '''
    # get the top b and top y scores from the database
    top_b, top_y = search.search_proteins(spectrum, database, TOP_N)
    # try and find some alignment from these scores
    alignments = aligners.align_spectrum_by_protein_ions(spectrum['spectrum'], top_b, top_y)
    # TODO: if alignments don't return anything promising, attempt some sort of hybrid alignment
    # for now, just take the top by and the top y and return them as alignment 1 and 2
    alignments = [x for _, x in alignments.items()]
    if len(alignments) < n:
        
        num_ions = ceil((n-len(alignments))/2)
        b_side = [x for _, x in top_b.items()][:num_ions]
        y_side = [x for _, x in top_y.items()][:num_ions]
        alignments += b_side + y_side
        
    return {'spectrum': spectrum['spectrum'], 'scan_no': spectrum['scan_no'], 'alignments': alignments}


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
            entry = id_spectrum(spectrum, database)
            entry_name = '{}_{}'.format(spectrum_file, spectrum['scan_no'])
            results[entry_name] = entry

    return results