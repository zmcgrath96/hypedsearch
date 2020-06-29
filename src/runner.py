'''
exec.py

Author: Zachary McGrath 
Date: 6 April 2020

Executor for the program
In charge of the flow of the program
'''
from os import walk
from src.identfication import id_spectra
from src import summary

def run(args: dict) -> None:
    '''
    Executing function for the program

    Inputs:
        args:   object arguments from main. Should be validated in main. Attributes of args:
            spectra_folder:     (str) full path the the directory containing all spectra files
            database_file:      (str) full path to the .fasta database file
            output_dir:         (str) full path the the directory to save output to
            min_peptide_len:    (int) minimum peptide length to consider
            max_peptide_len:    (int) maximum peptide length to consider
            tolerance:          (float) the ppm tolerance to allow in search
            verbose:            (bool) extra printing
            scoring_alg:        (str) scoring algorithm to use
            DEBUG:              (bool) debuging print messages. Default=False
    Outputs:
        None
    '''
    # get all the spectra file names
    spectra_files = []
    for (_, _, filenames) in walk(args['spectra_folder']):
        for fname in filenames:
            if '.mzml' not in fname.lower():
                continue
            spectra_files.append(args['spectra_folder'] + fname)
        break

    matched_spectra = id_spectra.id_spectra(
        spectra_files, args['database_file'], 
        min_peptide_len=args['min_peptide_len'], 
        max_peptide_len=args['max_peptide_len'], 
        ppm_tolerance=args['tolerance'], 
        verbose=True, 
        scoring_alg=args['scoring_alg'], 
        DEBUG=args['DEBUG']
    )
    print('\nFinished search. Writting results to {}...'.format(args['output_dir']))
    summary.generate(matched_spectra, args['output_dir'])
    
    