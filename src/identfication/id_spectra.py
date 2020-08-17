from src.file_io import mzML
from src.objects import Database, Spectrum, Alignments
from src.sequence.gen_spectra import gen_spectrum, gen_min_ordering
from src.identfication.alignment import attempt_alignment
from src.utils import ppm_to_da
from src.scoring import scoring
from src import database

from math import ceil

import bisect

def reduce_search_space(db: Database, reduced_spectra: list, ppm_tol: int) -> None:
    '''
    Look through the list of all masses in the experiment and find the interesting masses
    to keep in a graph for quicker searches later

    Inputs:
        db:                 (Database) the namedtuple instance that holds the tree, graphs, and proteins
        reduced_spectra:    (list) the list of all spectra to search the database for. MUST BE SORTED
        ppm_tol:            (int) the tolerance in parts per million to allow in searching for a mass
    Outputs:
        None
    '''
    # easiest way to keep track of unique kmers for each ion type
    b_hits = {}
    y_hits = {}

    # the maximum lenth peptide we could possibly find. Calculation is:
    # the ceiling of ( the maximum observed mass / the smallest singly charged mass (G))
    max_len = int(ceil(reduced_spectra[-1]/ 57.021464))

    def find_kmer_hits(kmer: str, prot_name: str) -> None:
        '''
        Internal function for adding both ions and charges to the tree 
        and to the dictionaries too keep track of sequences
        '''
        for ion in 'by':
   
            spec = gen_spectrum(kmer, ion=ion)['spectrum']

            for c, mass in enumerate(spec):

                # calculate the upper and lower bounds in the search
                da_tol = ppm_to_da(mass, ppm_tol)
                lb = mass - da_tol
                ub = mass + da_tol

                # get the index of the value to the left of any existing values of 
                # the lower bound. Any value after this postion is greater or equal to lower bound
                beginning_entry = bisect.bisect_left(reduced_spectra, lb)

                # see if the NEXT value is in the range. If so, keep the kmer
                if beginning_entry + 1 < len(reduced_spectra) and reduced_spectra[beginning_entry] <= ub:
                    # add the correct kmer to the correct dictionary and to the tree
                    if ion == 'b':
                        b_hits[kmer[:c+1]] = None
                        db.tree.insert(prot_name, kmer[:c+1])
                    else:
                        y_hits[kmer[-c-1:]] = None
                        db.tree.insert(prot_name, kmer[-c-1:])

    # instead of calling len all the time, call it once
    plen = len(db.proteins)

    print('Indexing potential hits...')

    # go through each protein and 
    for i, (prot_name, prot_entry) in enumerate(db.proteins.items()):
        
        print(f'\rOn protein {i+1}/{plen} [{int((i+1) * 100 / plen)}%]', end='')
        
        # iteratively go from 0, 
        for j in range(1, max_len):
            kmer = prot_entry.sequence[:j]
            find_kmer_hits(kmer, prot_name)
        
        for j in range(len(prot_entry.sequence) - max_len):
            kmer = prot_entry.sequence[j:j+max_len]
            find_kmer_hits(kmer, prot_name)     
            
        for j in range(len(prot_entry.sequence) - max_len, len(prot_entry.sequence)):
            kmer = prot_entry.sequence[j:]
            find_kmer_hits(kmer, prot_name)

    # go through all the b_hits and y_hits and add them to the corresponding graph
    print()
    blen = len(b_hits)
    for i, kmer in enumerate(sorted(b_hits.keys(), key=gen_min_ordering)):
        print(f'\rIndexing b k-mer {i+1}/{blen} [{int(100 * (i+1)/blen)}%]', end='')
        bs = gen_spectrum(kmer, ion='b', charge=1)['spectrum']
        bd = gen_spectrum(kmer, ion='b', charge=2)['spectrum']
        db.b_dawg.insert(bs, bd, kmer)
    db.b_dawg.finish()

    print()
    ylen = len(y_hits)
    for i, kmer in enumerate(sorted(y_hits.keys(), key=lambda x: gen_min_ordering(x[::-1]))):
        print(f'\rIndexing y k-mer {i+1}/{ylen} [{int(100 * (i+1)/ylen)}%]', end='')
        ys = gen_spectrum(kmer, ion='y', charge=1)['spectrum']
        yd = gen_spectrum(kmer, ion='y', charge=2)['spectrum']
        db.y_dawg.insert(ys, yd, kmer[::-1])
    db.y_dawg.finish()

    print('\nDone.')


def load_spectra(spectra_files: list, peak_filter=0, relative_abundance_filter=0) -> (list, list):
    '''
    Load all the spectra files into memory and merge all 
    spectra into one massive list for reduction of the search space

    Inputs:
        spectra_files:  (list) full paths to all spectra files
    kwargs:
        peak_filter:    (int) the top X most abundant spectra to keep. Default=0
        relative_abundance_filter:  (float) the percentage (0, 1) of the total abundance 
                                    a peak must have to be considered significant. Default=0
    Outputs:
        (list, list) all spectra in memory of Spectrum namedtuples, a list of floats of all masses
    '''
    # the single dimension of all the rounded spectra to use for reducing the search space
    linear_spectra = []

    # a list of all Spectrum namedtuples
    all_spectra = []

    # go through each spectra file and load them into memory
    for spectra_file in spectra_files:
        these_spectra = mzML.read(spectra_file, peak_filter=peak_filter, relative_abundance_filter=relative_abundance_filter)
        all_spectra += these_spectra

        # go through each mass of each spectrum, load it into memory, 
        # round the numbers to 3 decimal places for easier, and append to linear_spectra
        linear_spectra += list(set([
            round(x, 3) for spectrum in these_spectra for x in spectrum.spectrum
        ]))

    return (all_spectra, sorted(list(set(linear_spectra))))


def id_spectrum(db: Database, spectrum: Spectrum, ppm_tol: int, num_alignments=3) -> Alignments:
    '''
    Make num_alignments SequenceAlignments or HybridSequenceAlignments for a given spectrum

    Inputs:
        db:             (Database) the database to search
        spectrum:       (Spectrum) the spectrum namedtuple to make an alignment for
        ppm_tol:        (int) the tolerance in parts per million to allow when trying to match a mass
    kwargs:
        num_alignemnts: (int) the number of alignemnts to create for the spectrum
    Outputs:
        (Alignments) the namedtuple holding the spectrum and all of the attempted alignments
    '''
    # get the b and y hits for the spectrum
    b_kmers = db.b_dawg.fuzzy_search(spectrum.spectrum, 1, ppm_tol)
    sorted_b_results = sorted(
        [(kmer, scoring.score_subsequence(spectrum.spectrum, kmer, ppm_tol)[0]) for kmer in b_kmers],
         key=lambda x: x[1], reverse=True
    )
    max_score = sorted_b_results[0][1]
    filtered_b_results = [x[0] for x in sorted_b_results if x[1] == max_score]

    y_kmers = db.y_dawg.fuzzy_search(spectrum.spectrum, 1, ppm_tol)
    sorted_y_results = sorted(
        [(kmer, scoring.score_subsequence(spectrum.spectrum, kmer, ppm_tol)[1]) for kmer in y_kmers], 
        key=lambda x: x[1], 
        reverse=True
    )
    max_score = sorted_y_results[0][1]
    filtered_y_results = [x[0] for x in sorted_y_results if x[1] == max_score]

    return attempt_alignment(
        spectrum,
        db,
        filtered_b_results,
        filtered_y_results,
        ppm_tolerance=ppm_tol,
        n=num_alignments
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
    # build/load the database
    verbose and print('Loading database...')
    db: Database
    db = database.build(database_file)
    verbose and print('Done')
    
    # load all of the spectra
    verbose and print('Loading spectra...')
    spectra, reduced_spectra = load_spectra(
        spectra_files, 
        peak_filter=peak_filter, 
        relative_abundance_filter=relative_abundance_filter
    )
    verbose and print('Done')

    # build tree and graphs
    reduce_search_space(db, reduced_spectra, ppm_tolerance)

    # create an alignment for each spectrum
    results = {}
    spectrum: Spectrum

    num_spec = len(spectra)
    for i, spectrum in enumerate(spectra):
        print(f'On spectrum {i+1}/{num_spec} [{int(100 * (i+1) / num_spec)}%]')
        spectrum_key = f'{spectrum.file_name}_{spectrum.scan_number}'
        results[spectrum_key] = id_spectrum(db, spectrum, ppm_tolerance)

    return results