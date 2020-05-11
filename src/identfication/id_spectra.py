from src.file_io import mzML
from src.alignment import search, aligners
from src.database.database import Database
from src.modules.spectrum import Spectrum
from src.modules.alignment import Aligner
from math import ceil

############## Constants ##############
TOP_N = 5
#######################################

def id_spectrum(spectrum: Spectrum, database: Database, n=TOP_N, min_peptide_len=5, max_peptide_len=20) -> dict:
    '''
    Run a scoring and alignment on a specific spectrum 

    Inputs:
        spectrum:   Spectrum object
        database:   Database class instance
    kwargs:
        n:          number of alignments to return for each spectrum. Default=3
        min_peptide_len: minimum peptide length to consider. Default=5
        max_peptide_len: maximum peptide length to consider. Default=20
    Outputs:
        dict of the following form
        {
            'spectrum': list,
            'alignments': list of dicts,
            'scan_no': int
        }
    '''
    # get the top b and top y scores from the database
    top_b, top_y = search.search_database(spectrum, database, TOP_N)
    # top b and y are dictionaries with rank as key and values are Scores objects
    # try and find some alignment from these scores
    aligner = Aligner(spectrum)
    aligner.add_scores([x for _, x in top_b.items()], 'b')
    aligner.add_scores([x for _, x in top_y.items()], 'y')
    aligner.make_alignment()
        
    return aligner


def id_spectra(spectra_files: list, database_file: str, verbose=True, min_peptide_len=5, max_peptide_len=20) -> dict:
    '''
    Run a scoring and alignment on each spectra passed in

    Inputs:
        spectra_files:  list of strings of file names of spectra
        database_file:  string full path to a .fasta database
    kwargs: 
        verbose:        bool whether or not to print messages. Default=True
    Outputs:
        dict containing the results. All information is keyed by the spectrum file name with scan number appended and the values are Aligner objects
    '''
    verbose and print('Loading database...')
    database = Database(database_file, is_uniprot=True, kmer_size=min_peptide_len)
    verbose and print('Done. Indexing database...')
    database.index()
    verbose and print('Done.')
    verbose and print('Number of {}-mers found in the database: {}'.format(min_peptide_len, len(database.metadata.keys())))
    results = {}
    # go through all of the mzml files
    for i, spectrum_file in enumerate(spectra_files):
        print('Analyzing spectra file {}/{}[{}%]\n'.format(i, len(spectra_files), int(float(i)/float(len(spectra_files)) * 100)))

        spectra = mzML.read(spectrum_file)
        # go through each spectrum in the mzml file
        for j, spec in enumerate(spectra):
            print('Analyzing spectrum {}/{}[{}%]\r'.format(j, len(spectra), int(float(j)/float(len(spectra)) * 100)), end='')
            # make a spectrum object 
            spectrum = Spectrum(spec['spectrum'], spec['scan_no'], spec['level'], spec['precursor'], spec['abundance'])
            entry = id_spectrum(spectrum, database, min_peptide_len=min_peptide_len, max_peptide_len=max_peptide_len)
            entry_name = '{}_{}'.format(spectrum_file, spectrum.scan_number)
            results[entry_name] = entry

    return results