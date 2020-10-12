from src.file_io import mzML
from src.objects import Database, Spectrum, Alignments
from src.cppModules import gen_spectra, map_spectra_masses
from src.alignment import alignment
from src.utils import ppm_to_da, to_percent, overlap_intervals, predicted_len
from src.scoring import scoring, mass_comparisons
from src import database

from math import ceil
from collections import defaultdict, namedtuple
from operator import itemgetter
from typing import Iterable

import bisect
import time
import array as arr

BasicProtein = namedtuple('BasicProtein', ['name', 'sequence'], defaults=['', ''])

# top results to keep for creating an alignment
TOP_X = 50
BATCH_SIZE = 280

def hashable_boundaries(boundaries: list) -> str:
    return '-'.join([str(x) for x in boundaries])

def load_spectra(spectra_files: list, ppm_tol: int, peak_filter=0, relative_abundance_filter=0) -> (list, list, dict):
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
        (list, dict, list) all spectra in memory of boundaries namedtuples, a list of floats of all masses
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

    # turn the all spectra list into a list of boundaries
    def make_boundaries(mz):
        da_tol = ppm_to_da(mz, ppm_tol)
        return [mz - da_tol, mz + da_tol]

    boundaries = [make_boundaries(mz) for mz in linear_spectra]

    # make overlapped boundaries larger boundaries
    boundaries = overlap_intervals(boundaries)

    # make a mapping for mz -> boundaries
    b_i, s_i = 0, 0
    mz_mapping = {}
    while s_i < len(linear_spectra):
        
        # if the current boundary encapsulates s_i, add to list
        if boundaries[b_i][0] <= linear_spectra[s_i] <= boundaries[b_i][1]:
            mz_mapping[linear_spectra[s_i]] = b_i 
            s_i += 1

        # else if the s_i < boundary, increment s_i
        elif linear_spectra[s_i] < boundaries[b_i][0]:
            s_i += 1

        # else if s_i > boundary, incrment b_i
        elif linear_spectra[s_i] > boundaries[b_i][1]:
            b_i += 1

    return (all_spectra, boundaries, mz_mapping)

def match_masses(spectra_boundaries: list, db: Database) -> (dict, dict, Database):
    '''
    Take in a list of boundaries from observed spectra and return a b and y
    dictionary that maps boundaries -> kmers

    Inputs:
        spectra_boundaries: (list) lists with [lower_bound, upper_bound]
        db:                 (Database) entries to look for kmers in
    Outputs:
        (dict, dict) mappings from boundaries -> kmers for (b, y) ions
    '''

    '''
    We need to create mappings from boundaries -> kmers for both b and y ions
    Emprically we found that ~300 proteins worth of kmers and floats takes ~4GB RAM
    So in order to keep in down we need to batch it. Heres what we'll do
    
      1. Define a constant above for the number of proteins to have per batch
      2. For the number of proteins in a batch
        1. Pass those proteins to make_database_set (returns a list of b, y kmers with masses and indices)
        2. Pass those lists ot merge
        3. The results of that are mappings of boundaries -> kmers, add those to some constant mapping we have
        4. Delete those old lists to clear memory
        5. Repeat
    
    Adding each batch's results to our mapping will create a large mapping. We will return
    that mapping and the updated database

    '''
    # keep track of all of the good mass matches and kmers
    matched_masses_b, matched_masses_y, kmer_set = defaultdict(list), defaultdict(list), defaultdict(list)

    # estimate the max len
    max_len = predicted_len(spectra_boundaries[-1][1])

    # calc the number of batches needed
    num_batches = ceil(len(db.proteins) / BATCH_SIZE)

    # create batches of proteins in the form of (prot name, prot entry)
    kv_prots = [BasicProtein(k, v.sequence) for k, v in db.proteins.items()]
    batched_prots = [kv_prots[i*BATCH_SIZE:(i+1)*BATCH_SIZE] for i in range(num_batches)]

    # go through each batch set and create the list representation, merge, and keep good prots
    for batch_num, batch_set in enumerate(batched_prots):

        print(f'On batch {batch_num + 1}/{num_batches}\n', end='')
        batch_start_time = time.time()

        # create our list representation
        matched_masses_b_batch, matched_masses_y_batch, batch_kmer_set = map_spectra_masses.map_boundaries(
            spectra_boundaries, 
            batch_set, 
            max_len
        )

        # add these these hits to our function scoped set of masses and to the kmer set
        for k, v in matched_masses_b_batch.items():
            matched_masses_b[k] += v 

        for k, v in matched_masses_y_batch.items():
            matched_masses_y[k] += v 

        for k, v in batch_kmer_set.items():
            kmer_set[k] += v

        print(f'\nBatch time is : {time.time() - batch_start_time} seconds')

    #update kmers in db to the kmer set
    db = db._replace(kmers=kmer_set)

    return (matched_masses_b, matched_masses_y, db)


def id_spectrum(
    spectrum: Spectrum, 
    boundaries: list, 
    db: Database,
    mz_mapping: dict, 
    matched_masses_b: dict, 
    matched_masses_y: dict,
    ppm_tolerance: int
    ) -> Alignments:
    '''
    Create an alignment for a spectrum

    Inputs:
        spectrum:           (Spectrum) the spectrum to create an alignment for
        boundaries:         (list) the list of boundaries for all observed spectra
        db:                 (Database) holding all the records
        mz_mapping:         (dict) the mapping from mass -> boundaries needed to retrieve kmers for a mass (boundary)
        matched_masses_b:   (dict) the mapping from boundaries -> kmers for b ions
        matched_masses_y:   (dict) the mapping from boundaries -> kmers for y ions
        ppm_tolerance:      (int) the tolerance to allow for scoring algs
    Outputs:
        (Alignments) the created alignments for this spectrum
    '''
    b_hits, y_hits = [], []
    for mz in spectrum.spectrum:

        # get the correct boundary
        mapped = mz_mapping[mz]
        b = boundaries[mapped]
        b = hashable_boundaries(b)

        if b in matched_masses_b:
            b_hits += matched_masses_b[b]

        if b in matched_masses_y:
            y_hits += matched_masses_y[b]

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
    db: Database
    db = database.build(database_file)
    verbose and print('Done')
    
    # load all of the spectra
    verbose and print('Loading spectra...')
    spectra, boundaries, mz_mapping = load_spectra(
        spectra_files, 
        ppm_tolerance,
        peak_filter=peak_filter, 
        relative_abundance_filter=relative_abundance_filter
    )
    verbose and print('Done')

    # get the boundary -> kmer mappings for b and y ions
    print('Indexing potential hits...')
    matched_masses_b, matched_masses_y, db = match_masses(boundaries, db)
    print('Done')

    # keep track of the alingment made for every spectrum
    results = {}
    for i, spectrum in enumerate(spectra):
        print(f'\rCreating an alignment for {i+1}/{len(spectra)} [{to_percent(i, len(spectra))}%]', end='')
        results[i] = id_spectrum(spectrum, boundaries, db, mz_mapping, matched_masses_b, matched_masses_y, ppm_tolerance)
        
    return results