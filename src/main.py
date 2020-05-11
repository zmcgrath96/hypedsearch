'''
main.py

Author: Zach McGrath
Date: 6 April 2020

Description:
Main file for hypedsearch. Handles input parameters and flow of the program
'''
import argparse
from src.utils import utils
import sys
from src import runner

##############################################################

def main(args: object) -> None:
    ############## Argument checking ################
    if not utils.is_dir(args.spectra_folder):
        print('Error: {} is not a real path. Path to directory with spectra files is necessary.')
        sys.exit(0)
    if not utils.is_fasta(args.database_file) or not utils.is_file(args.database_file):
        print('Error: {} is not a valid .fasta file. .fasta file needed.')
        sys.exit(0)

    output_dir = utils.make_valid_dir_string(args.save_dir)
    utils.make_dir(output_dir)
    ###################  Run  #######################
    arguments = {
        'spectra_folder': args.spectra_folder,
        'database_file': args.database_file,
        'output_dir': output_dir,
        'min_peptide_len': args.min_peptide_len,
        'max_peptide_len': args.max_peptide_len,
    }
    runner.run(arguments)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Entry file for the database experiments')

    parser.add_argument('spectra_folder', type=str, metavar='SF', help='Path to folder containing spectra files.')
    parser.add_argument('database_file', type=str, metavar='DB', help='Path to .fasta file containing proteins')
    parser.add_argument('--output-dir', dest='save_dir', type=str, default='~/', help='Directory to save all figures. Default=~/')
    parser.add_argument('--min-peptide-len', dest='min_peptide_len', type=int, default=5, help='Minimum peptide length to consider. Default=5')
    parser.add_argument('--max-peptide-len', dest='max_peptide_len', type=int, default=20, help='Maximum peptide length to consider. Default=20')
    args = parser.parse_args()
    main(args)