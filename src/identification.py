from src.file_io import mzML
from src.objects import Spectrum, Alignments
from src.cppModules import gen_spectra
from src.alignment import alignment
from src.utils import to_percent 
from src.scoring import scoring, mass_comparisons
from src.database import Database

# top results to keep for creating an alignment
TOP_X = 50

def load_spectra(spectra_files: list, peak_filter=0, relative_abundance_filter=0) -> (list, dict):
    '''
    Load all the spectra files into memory and merge all 
    spectra into one massive list for reduction of the search space

    indicesnputs:
        spectra_files:  (list) full paths to all spectra files
    kwargs:
        peak_filter:    (int) the top X most abundant spectra to keep. Default=0
        relative_abundance_filter:  (float) the percentage (0, 1) of the total abundance 
                                    a peak must have to be considered significant. Default=0
    Outputs:
        (list, list) list of Spectrum objects, list of all masses 
    '''
    # the single dimension of all the rounded spectra to use for reducing the search space
    linear_spectra = []

    # a list of all boundaries namedtuples
    all_spectra = []

    # go through each spectra file and load them into memory
    for spectra_file in spectra_files:
        these_spectra = mzML.read(spectra_file, peak_filter=peak_filter, relative_abundance_filter=relative_abundance_filter)
        all_spectra += these_spectra

        # go through each mass of each s, load it into memory, 
        # round the numbers to 3 decimal places for easier, and append to linear_spectra
        linear_spectra += list(set([
            x for spectrum in these_spectra for x in spectrum.spectrum
        ]))

    # sort the linear spectra
    linear_spectra.sort()

    return (all_spectra, linear_spectra)

def id_spectrum(spectrum: Spectrum, db: Database, ppm_tolerance: int) -> Alignments:
    '''
    Create an alignment for a spectrum

    Inputs:
        spectrum:           (Spectrum) the spectrum to create an alignment for
        db:                 (Database) holding all the records
        ppm_tolerance:      (int) the tolerance to allow for scoring algs
    Outputs:
        (Alignments) the created alignments for this spectrum
    '''
    b_hits, y_hits = [], []
    for mz in spectrum.spectrum:

        mz_b_hits, mz_y_hits = db.get_mass_tags(mz)
        b_hits += mz_b_hits
        y_hits += mz_y_hits

    # remove any duplicates
    b_hits = list(set(b_hits))
    y_hits = list(set(y_hits))

    # score and sort these results
    b_results = sorted([
        (
            kmer, 
            mass_comparisons.optimized_compare_masses(spectrum.spectrum, gen_spectra.gen_spectrum(kmer, ion='b'))
        ) for kmer in b_hits], 
        key=lambda x: (x[1], 1/len(x[0])), 
        reverse=True
    )
    y_results = sorted([
        (
            kmer, 
            mass_comparisons.optimized_compare_masses(spectrum.spectrum, gen_spectra.gen_spectrum(kmer, ion='y'))
        ) for kmer in y_hits], 
        key=lambda x: (x[1], 1/len(x[0])), 
        reverse=True
    )

    # filter out the results
    # 1. take all non-zero values 
    # 2. either take the TOP_X or if > TOP_X have the same score, all of those values
    filtered_b, filtered_y = [], []

    # find the highest b and y scores
    max_b_score = max([x[1] for x in b_results])
    max_y_score = max([x[1] for x in y_results])

    # count the number fo kmers that have the highest value
    num_max_b = sum([1 for x in b_results if x[1] == max_b_score])
    num_max_y = sum([1 for x in y_results if x[1] == max_y_score])

    # if we have more than TOP_X number of the highest score, take all of them
    keep_b_count = max(TOP_X, num_max_b)
    keep_y_count = max(TOP_X, num_max_y)

    # take the afformentioned number of results that > than zero
    filtered_b = [x[0] for x in b_results[:keep_b_count] if x[1] > 0]
    filtered_y = [x[0] for x in y_results[:keep_y_count] if x[1] > 0]

    # create an alignment for the spectrum
    return alignment.attempt_alignment(
        spectrum, 
        db, 
        filtered_b, 
        filtered_y, 
        ppm_tolerance=ppm_tolerance, 
        n=5
    )


def id_spectra(
    spectra_files: list, 
    database_file: str, 
    verbose=True, 
    min_peptide_len=5, 
    max_peptide_len=20, 
    result_count=3, 
    peak_filter=0, 
    relative_abundance_filter=0.0,
    ppm_tolerance=20, 
    precursor_tolerance=1, 
    DEBUG=False
) -> dict:
    '''
    Run a scoring and alignment on each spectra passed in and give spectra a sequence of 
    Amino Acids that best describes the spectrum

    indicesnputs:
        spectra_files:          (list of strings) of file names of spectra
        database_file:          (string) full path to a .fasta database
    kwargs: 
        verbose:                (bool) whether or not to print messages. Default=True
        min_peptide_len:        (int) minimum length sequence to consider for alignment. Default=5
        max_peptide_len:        (int) maximum length sequence to consider for alignemtn. Default=20
        ppm_tolerance:          (int) tolerance for ppm to include in search. Default=20
        precursor_tolerance:    (float) the tolerance to allow when matching precusor masses. Default=1
    Outputs:
        dict containing the results. 
        All information is keyed by the spectrum file name with scan number appended 
        and the values are list of boundariesequenceAligment objects
    '''
    # build/load the database
    verbose and print('Loading database...')
    db = Database(database_file)
    verbose and print('Done')

    print(f'Number of proteins: {len(db.proteins)}')
    
    # load all of the spectra
    verbose and print('Loading spectra...')
    spectra, linear_spectra = load_spectra(
        spectra_files, 
        peak_filter=peak_filter, 
        relative_abundance_filter=relative_abundance_filter
    )
    verbose and print('Done')

    # index the database
    db.reduce_search_space(linear_spectra, ppm_tolerance)

    # keep track of the alingment made for every spectrum
    results = {}
    for i, spectrum in enumerate(spectra):
        print(f'\rCreating an alignment for {i+1}/{len(spectra)} [{to_percent(i, len(spectra))}%]', end='')
        results[i] = id_spectrum(spectrum, db, ppm_tolerance)
        
    return results