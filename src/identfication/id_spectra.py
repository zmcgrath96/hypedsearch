from src.file_io import mzML
from src.identfication.alignment import attempt_alignment
from src import database
from src.objects import Spectrum, MassSequence, KmerMassesResults, Alignments, Database
from src.scoring import mass_comparisons
from src.sequence.gen_spectra import gen_spectrum

from collections import defaultdict
import math 
from bisect import bisect


def id_spectrum(
    spectrum: Spectrum, 
    db: Database, 
    min_peptide_len: int, 
    n=3, 
    ppm_tolerance=20, 
    precursor_tolerance=1,
    scoring_alg='ibb', 
    DEBUG=False
) -> Alignments:
    '''
    Run an alignemnt in the form of an Amino Acid sequence with a score to the caller for this spectrum.
    The top n results are returned. If no alignments can be created, None is returned.

    Inputs:
        spectrum:               (Spectrum) the spectrum to id
        kmermasses:             (KmerMasses) hash tables with al the necessary dictionaries
    kwargs:
        n:                      (int) number of alignments to return for each spectrum. Default=3
        ppm_tolerance:          (int) the ppm tolerance to allow when searching. Default=20
        precursor_tolerance:    (float) the tolerance to allow when matching precusor masses. Default=1
        scoring_alg:            (str) the name of the scoring algoirhtm to use. Options are 'bb', 'ion', 'ibb'. Default=bb
    Outputs:
        Alignments namedtuple or None
    '''
    # search the mass tables
    bs, bd, ys, yd = database.search(db, spectrum, 'bs', 20), \
                    database.search(db, spectrum, 'bd', 20), \
                    database.search(db, spectrum, 'ys', 20), \
                    database.search(db, spectrum, 'yd', 20)
    
    # put the results into a structrue
    hits = KmerMassesResults(bs, bd, ys, yd)

    # if we get no hits whatsoever, return None
    if all([len(x) == 0 for x in hits]):
        return None
        
    # attempt alignments
    a = attempt_alignment(
        spectrum, 
        db, 
        hits, 
        min_peptide_len, 
        ppm_tolerance=ppm_tolerance, 
        scoring_alg=scoring_alg, 
        precursor_tolerance=precursor_tolerance, 
        DEBUG=DEBUG
    )
    
    # return the alignments in a structure
    return Alignments(spectrum, a)

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
    scoring_alg='ibb', 
    DEBUG=False
) -> dict:
    '''
    Run a scoring and alignment on each spectra passed in and give spectra a sequence of 
    Amino Acids that best describes the spectrum

    Inputs:
        spectra_files:          (list of strings) of file names of spectra
        database_file:          (string) full path to a .fasta database
    kwargs: 
        verbose:                (bool) whether or not to print messages. Default=True
        min_peptide_len:        (int) minimum length sequence to consider for alignment. Default=5
        max_peptide_len:        (int) maximum length sequence to consider for alignemtn. Default=20
        ppm_tolerance:          (int) tolerance for ppm to include in search. Default=20
        precursor_tolerance:    (float) the tolerance to allow when matching precusor masses. Default=1
        scoring_alg:            (str) the scoring algorithm to use. Either 'bb' or 'ion'. Default=bb
    Outputs:
        dict containing the results. 
        All information is keyed by the spectrum file name with scan number appended 
        and the values are list of SequenceAligment objects
    '''

    # load the datbase into memory
    verbose and print('Loading database...')
    db = database.build_or_load_db(Database(database_file, min_peptide_len, max_peptide_len, verbose))
    verbose and print('Done.')

    # keep track of the results
    results = {}

    # go through all of the spectra files
    for i, spectrum_file in enumerate(spectra_files):

        verbose and print('Analyzing spectra file {}/{}[{}%]\n'.format(i + 1, len(spectra_files), int(float(i)/float(len(spectra_files)) * 100)))

        # load the spectra into memory
        spectra = mzML.read(spectrum_file, peak_filter=peak_filter, relative_abundance_filter=relative_abundance_filter)

        # go through each spectrum in the mzml file
        for j, spec in enumerate(spectra):

            verbose and print('Analyzing spectrum {}/{}[{}%]\r'.format(j + 1, len(spectra), int(float(j)/float(len(spectra)) * 100)), end='')            
            
            # align this spectrum
            aligned_spectrum = id_spectrum(
                spec,
                 db,
                 min_peptide_len,
                 result_count,
                 ppm_tolerance=ppm_tolerance,
                 precursor_tolerance=precursor_tolerance,
                 scoring_alg=scoring_alg,
                DEBUG=DEBUG
            )
            
            # if an alignment cannot be created, continue
            if aligned_spectrum is None:
                continue

            # save the results in the dictionary for now
            entry_name = '{}_{}'.format(spectrum_file, spec.scan_number)
            results[entry_name] = aligned_spectrum

    return results