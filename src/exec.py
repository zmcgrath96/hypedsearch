'''
exec.py

Author: Zachary McGrath 
Date: 6 April 2020

Executor for the program
In charge of the flow of the program
'''
from os import walk
from src.identfication import id_spectra

def execute(args: dict) -> None:
    '''
    Executing function for the program

    Inputs:
        args:   object arguments from main. Should be validated in main. Attributes of args:
            'spectra_folder':   string full path the the directory containing all spectra files
            'database_file':    string full path to the .fasta database file
            'output_dir':       string full path the the directory to save output to
    Outputs:
        None
    '''
    # get all the spectra file names
    spectra_files = []
    for (_, _, filenames) in walk(args['spectra_folder']):
        for fname in filenames:
            spectra_files.append(args['spectra_folder'] + fname)
        break

    matched_spectra = id_spectra.id_spectra(spectra_files, args['database_file'])
    
    