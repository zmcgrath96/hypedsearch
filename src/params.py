########################## PROGRAM IO ##########################
# full path to the spectra folder containing all spectra files
#SPECTRA_FOLDER = '/Users/zacharymcgrath/Desktop/nod2 data/filteredSpec/'
SPECTRA_FOLDER = '/Users/zacharymcgrath/Desktop/nod2 data/single/'
# full path to the .fasta database file
#DATABASE_FILE = '/Users/zacharymcgrath/Desktop/nod2 data/filteredNOD2.fasta'
DATABASE_FILE = '/Users/zacharymcgrath/Desktop/raw inputs/NOD2_E3/filtered_mouse_database.fasta'
# full path to the output directory
#OUTPUT_DIRECTORY = '/Users/zacharymcgrath/Desktop/Experiment output/filtered_NOD2_E3_SEP_22/'
OUTPUT_DIRECTORY = '/Users/zacharymcgrath/Desktop/raw inputs/NOD2_E3/output/'

########################## SEARCH PARAMETERS ##########################
# minimum length peptide to look for
MIN_PEPTIDE_LEN = 3
# maximum length peptide to look for
MAX_PEPTIDE_LEN = 30
# tolerance (in ppm) to allow when matching m/z peaks
PPM_TOLERANCE = 20
# tolerance (in Da) to allow when matching precursor masses
PRECURSOR_TOLERANCE = 10
#-----------------------------------------------------------
# NOTE: when using the filtering options below (either the peak
# filtering or the abundance filtering), only set 1 of the values.
# If both are set, then it defaults to the value set for peak_filter
#-----------------------------------------------------------
# the number of most abundant (intense) peaks to allow when filtering
# input spectra
PEAK_FILTER = 25
# the minimum allowed percentage of the total intensity a peak 
# is allowed to have to pass the filter. Values should be 
# in the range (0, 1) 
RELATIVE_ABUNDANCE_FILTER = 0

########################## EXPERIMENT PARAMETERS ##########################
# in vitro digest performed on the sample. Leave '' if none
DIGEST = 'trypsin'

########################## RUN PARAMETERS ##########################
# printing to the console during the run
VERBOSE = True
# Extra extra printing options for debugging the application
DEBUG = True
# the number of cores to allow in the search. Should be a number 
CORES = 1
# the number of alignments to keep per spectrum
N = 5

######################### DEV PARAMS #########################
# "truth" set json file that has the "truth" sequence for every 
# spectrum. Should be keyed by some id param that is extracted by the tool 
# to check where the correct value falls off at. Format of the json is simple:
# {spectrum_id: {'sequence': str, 'hybrid': bool, 'parent': str}}
# if the param is left blank or the file cannot be found, it is not used
TRUTH_SET = '' #'/Users/zacharymcgrath/Desktop/Experiment output/fall_off/specmil_truth_set.json'
