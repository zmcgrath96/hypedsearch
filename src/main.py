'''
main.py

Author: Zach McGrath
Date: 6 April 2020

Description:
Main file for hypedsearch. Handles input parameters and flow of the program
'''
import argparse
from src import utils, runner
import sys

def stringtobool(s: str) -> bool:
    s = str(s)
    if s.lower() == 'false' or 'f' in s.lower():
        return False
    return True


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
        'tolerance': args.tolerance,
        'precursor_tolerance': args.precursor_tolerance,
        'verbose': stringtobool(args.verbose), 
        'scoring_alg': args.scoring_alg, 
        'peak_filter': args.peak_filter, 
        'relative_abundance_filter': args.rel_abund_filter,
        'DEBUG': False
    }
    runner.run(arguments)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Tool for identifying proteins, both hybrid and non hybrid from MS/MS data')

    parser.add_argument('spectra_folder', type=str, metavar='SF', help='Path to folder containing spectra files.')
    parser.add_argument('database_file', type=str, metavar='DB', help='Path to .fasta file containing proteins')
    parser.add_argument('--output-dir', dest='save_dir', type=str, default='~/', help='Directory to save all figures. Default=~/')
    parser.add_argument('--min-peptide-len', dest='min_peptide_len', type=int, default=5, help='Minimum peptide length to consider. Default=5')
    parser.add_argument('--max-peptide-len', dest='max_peptide_len', type=int, default=20, help='Maximum peptide length to consider. Default=20')
    parser.add_argument('--tolerance', dest='tolerance', type=int, default=20, help='ppm tolerance to allow in search. Deafult=20')
    parser.add_argument('--precursor-tolerance', dest='precursor_tolerance', type=float, default=1, help='The mass (in Da) tolerance to accept when matching precursor masses. Default=1')
    parser.add_argument('--peak-filter', dest='peak_filter', type=int, default=0, help='The number of peaks to take from a spectrum. The most abundant peaks will be taken. Leave blank if you want no filter or to use relative abundance filter. Defualt=0')
    parser.add_argument('--abundance_filter', dest='rel_abund_filter', type=float, default=0.0, help='Take only peaks from a spectrum where the abundance of the peak is >= the percentage give. Leave blank if you want no filter or to use peak filter. Default=0.0')
    parser.add_argument('--score', dest='scoring_alg', type=str, default='ibb', help='Scoring algorithm to use. Options are [bb, ion, ibb] for backbone, ion, and ion backbone respectively. Default=bb')
    parser.add_argument('--verbose', dest='verbose', type=bool, default=True, help='Extra printing to console during run. Default=True')
    args = parser.parse_args()
    main(args)