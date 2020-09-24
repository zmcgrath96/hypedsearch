'''
main.py

Author: Zach McGrath
Date: 6 April 2020

Description:
Main file for hypedsearch. Handles input parameters and flow of the program
'''
import argparse
from src import utils, runner, params
import sys

def stringtobool(s: str) -> bool:
    s = str(s)
    if s.lower() == 'false' or 'f' in s.lower():
        return False
    return True

def set_args(args) -> dict:

    # check if we use the params file or the user arguments
    use_params = stringtobool(args.params)

    spectra_folder = args.spectra_folder if not use_params else params.SPECTRA_FOLDER
    database_file = args.database_file if not use_params else params.DATABASE_FILE
    output_dir = args.output_dir if not use_params else params.OUTPUT_DIRECTORY
    min_peptide_len = args.min_peptide_len if not use_params else params.MIN_PEPTIDE_LEN
    max_peptide_len = args.max_peptide_len if not use_params else params.MAX_PEPTIDE_LEN
    ppm_tolerance = args.tolerance if not use_params else params.PPM_TOLERANCE
    precursor_tolerance = args.precursor_tolerance if not use_params else params.PRECURSOR_TOLERANCE
    verbose = stringtobool(args.verbose) if not use_params else params.VERBOSE
    peak_filter = args.peak_filter if not use_params else params.PEAK_FILTER
    relative_abundance_filter = args.rel_abund_filter if not use_params else params.RELATIVE_ABUNDANCE_FILTER
    debug = params.DEBUG

    ############## Argument checking ################
    if not utils.is_dir(spectra_folder):
        print(f'Error: {spectra_folder} is not a real path. Path to directory with spectra files is necessary.')
        sys.exit(0)
    if not utils.is_fasta(database_file) or not utils.is_file(database_file):
        print(f'Error: {database_file} is not a valid .fasta file. .fasta file needed.')
        sys.exit(0)

    # make the output directory
    output_dir = utils.make_valid_dir_string(output_dir)
    utils.make_dir(output_dir)

    return {
        'spectra_folder': spectra_folder,
        'database_file': database_file,
        'output_dir': output_dir,
        'min_peptide_len': min_peptide_len,
        'max_peptide_len': max_peptide_len,
        'tolerance': ppm_tolerance,
        'precursor_tolerance': precursor_tolerance,
        'verbose': verbose, 
        'peak_filter': peak_filter, 
        'relative_abundance_filter': relative_abundance_filter,
        'DEBUG': debug
    }

##############################################################

def main(args: object) -> None:
    # get the arguments 
    arguments = set_args(args)

    runner.run(arguments)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Tool for identifying proteins, both hybrid and non hybrid from MS/MS data')

    parser.add_argument('--spectra_folder', dest='spectra_folder', type=str, default='./', help='Path to folder containing spectra files.')
    parser.add_argument('--database_file', dest='database_file', type=str, default='./', help='Path to .fasta file containing proteins')
    parser.add_argument('--output-dir', dest='output_dir', type=str, default='~/', help='Directory to save all figures. Default=~/')
    parser.add_argument('--params', dest='params', type=bool, default=False, help='Use the params.py file adjacent to main.py instead of using command line arguments. Default=False')
    parser.add_argument('--min-peptide-len', dest='min_peptide_len', type=int, default=5, help='Minimum peptide length to consider. Default=5')
    parser.add_argument('--max-peptide-len', dest='max_peptide_len', type=int, default=20, help='Maximum peptide length to consider. Default=20')
    parser.add_argument('--tolerance', dest='tolerance', type=int, default=20, help='ppm tolerance to allow in search. Deafult=20')
    parser.add_argument('--precursor-tolerance', dest='precursor_tolerance', type=float, default=1, help='The mass (in Da) tolerance to accept when matching precursor masses. Default=1')
    parser.add_argument('--peak-filter', dest='peak_filter', type=int, default=0, help='The number of peaks to take from a spectrum. The most abundant peaks will be taken. Leave blank if you want no filter or to use relative abundance filter. Defualt=0')
    parser.add_argument('--abundance_filter', dest='rel_abund_filter', type=float, default=0.0, help='Take only peaks from a spectrum where the abundance of the peak is >= the percentage give. Leave blank if you want no filter or to use peak filter. Default=0.0')
    parser.add_argument('--verbose', dest='verbose', type=bool, default=True, help='Extra printing to console during run. Default=True')
    args = parser.parse_args()
    main(args)