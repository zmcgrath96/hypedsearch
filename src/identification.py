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
# BATCH_SIZE = 300

def hashable_boundaries(boundaries: list) -> str:
    return '-'.join([str(x) for x in boundaries])

# def merge(P: Iterable, indices: Iterable, kmers: Iterable, boundaries: Iterable):
#     b_i, p_i = 0, 0

#     matched_masses = defaultdict(list)

#     while b_i < len(boundaries) and p_i < len(P):

#         # if P[p_i] is in the boundary, keep track of it increment p_i
#         if boundaries[b_i][0] <= P[p_i] <= boundaries[b_i][1]:
#             matched_masses[hashable_boundaries(boundaries[b_i])] += kmers[indices[p_i - 1]:indices[p_i]]
#             p_i += 1

#         # if the upper boundary is less than the p_i, increment b_i
#         elif P[p_i] > boundaries[b_i][1]:
#             b_i += 1

#         # if the lower boundary is greater than the p_i, increment p_i
#         elif P[p_i] < boundaries[b_i][0]:
#             p_i += 1

#     return matched_masses

# def make_database_set(proteins: list, max_len: int) -> (list, list, list):
#     '''
#     Create parallel lists of (masses, index_maps, kmers) for the merge sort operation

#     Inputs:
#         proteins:   (list) protein entries of shape (name, entry) where entry has a .sequence attribute
#         max_len:    (int) the max length kmer to generate
#     Ouputs:
#         db_list_b:      (array) all b ion masses created in iteration through the proteins
#         index_list_b:   (array) the index mapping. Equal size to db_list_b where each entry maps 
#                                 idx(db_list_b) -> range of indices in kmer_list_b of kmers corresponding to masses
#         kmer_list_b:    (list) the list of kmers associated to db_list_b through index_list_b

#         db_list_y:      (array) all y ion masses created in iteration through the proteins
#         index_list_y:   (array) he index mapping. Equal size to db_list_y where each entry maps 
#                                 idx(db_list_y) -> range of indices in kmer_list_y of kmers corresponding to masses
#         kmer_list_y:    (list) the list of kmers associated to db_list_b through index_list_y

#         kmer_set:       (dict) mappings of kmers -> source proteins

#     '''
#     db_dict_b = defaultdict(set)
#     db_dict_y = defaultdict(set)
#     kmer_set = defaultdict(list)

#     # function to add all masses of b+, b++, y+, y++ ions to the correct dict for sorting
#     # and keep track of the source proteins
#     def add_all(kmer, prot_name):
#         for ion in 'by':
#             for charge in [1, 2]:
#                 spec = gen_spectra.gen_spectrum(kmer, ion=ion, charge=charge)

#                 for i, mz in enumerate(spec):
#                     kmer_to_add = kmer[:i+1] if ion == 'b' else kmer[-i-1:]
#                     r_d = db_dict_b if ion == 'b' else db_dict_y
#                     r_d[mz].add(kmer_to_add)
#                     kmer_set[kmer_to_add].append(prot_name)

#     plen = len(proteins)

#     # go through each protein and add all kmers to the correct dictionary for later sorting
#     for i, (prot_name, prot_entry) in enumerate(proteins):
        
#         print(f'\rOn protein {i+1}/{plen} [{int((i+1) * 100 / plen)}%]', end='')

#         for j in range(1, max_len):
#             kmer = prot_entry.sequence[:j]
#             add_all(kmer, prot_name)
            
#         for j in range(len(prot_entry.sequence) - max_len):
#             kmer = prot_entry.sequence[j:j+max_len]
#             add_all(kmer, prot_name)

#         for j in range(len(prot_entry.sequence) - max_len, len(prot_entry.sequence)):
#             kmer = prot_entry.sequence[j:]
#             add_all(kmer, prot_name)

#     print('\nSorting the set of protein masses...')
    
#     # arrays take less space than lists
#     db_list_b, index_list_b, kmer_list_b = arr.array('f'), arr.array('i'), []
#     db_list_y, index_list_y, kmer_list_y = arr.array('f'), arr.array('i'), []

#     # sort the keys to make sure that we are going through masses the correct way
#     sorted_keys = sorted(db_dict_b.keys())
#     for mz in sorted_keys:

#         # get all kmers associated with this mass, append it the b ion list, and keep track 
#         # of the kmers in the kmer list
#         kmers = db_dict_b[mz]
#         db_list_b.append(mz)
#         offset = 0 if not len(index_list_b) else index_list_b[-1]
#         index_list_b.append(len(kmers) + offset)
#         kmer_list_b += kmers

#     sorted_keys = sorted(db_dict_y.keys())
#     for mz in sorted_keys:

#         # get all kmers associated with this mass, append it the y ion list, and keep track 
#         # of the kmers in the kmer list
#         kmers = db_dict_y[mz]
#         db_list_y.append(mz)
#         offset = 0 if not len(index_list_y) else index_list_y[-1]
#         index_list_y.append(len(kmers) + offset)
#         kmer_list_y += kmers

#     return db_list_b, index_list_b, kmer_list_b, db_list_y, index_list_y, kmer_list_y, kmer_set

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

# def match_masses(spectra_boundaries: list, db: Database) -> (dict, dict, Database):
#     '''
#     Take in a list of boundaries from observed spectra and return a b and y
#     dictionary that maps boundaries -> kmers

#     Inputs:
#         spectra_boundaries: (list) lists with [lower_bound, upper_bound]
#         db:                 (Database) entries to look for kmers in
#     Outputs:
#         (dict, dict) mappings from boundaries -> kmers for (b, y) ions
#     '''

#     '''
#     We need to create mappings from boundaries -> kmers for both b and y ions
#     Emprically we found that ~300 proteins worth of kmers and floats takes ~4GB RAM
#     So in order to keep in down we need to batch it. Heres what we'll do
    
#       1. Define a constant above for the number of proteins to have per batch
#       2. For the number of proteins in a batch
#         1. Pass those proteins to make_database_set (returns a list of b, y kmers with masses and indices)
#         2. Pass those lists ot merge
#         3. The results of that are mappings of boundaries -> kmers, add those to some constant mapping we have
#         4. Delete those old lists to clear memory
#         5. Repeat
    
#     Adding each batch's results to our mapping will create a large mapping. We will return
#     that mapping and the updated database

#     '''
#     # keep track of all of the good mass matches and kmers
#     matched_masses_b, matched_masses_y, kmer_set = defaultdict(list), defaultdict(list), defaultdict(list)

#     # estimate the max len
#     max_len = predicted_len(spectra_boundaries[-1][1])

#     # calc the number of batches needed
#     num_batches = ceil(len(db.proteins) / BATCH_SIZE)

#     # create batches of proteins in the form of (prot name, prot entry)
#     kv_prots = [(k, v) for k, v in db.proteins.items()]
#     batched_prots = [kv_prots[i*BATCH_SIZE:(i+1)*BATCH_SIZE] for i in range(num_batches)]

#     # go through each batch set and create the list representation, merge, and keep good prots
#     for batch_num, batch_set in enumerate(batched_prots):

#         print(f'On batch {batch_num + 1}/{num_batches}\n', end='')

#         # create our list representation
#         batch_b_list, index_list_b, batch_kmer_b, batch_y_list, index_list_y, batch_kmer_y, batch_kmer_set = make_database_set(batch_set, max_len)

#         # find tha batch matched masses for both b and y ions
#         matched_masses_b_batch = merge(batch_b_list, index_list_b, batch_kmer_b, spectra_boundaries)
#         matched_masses_y_batch = merge(batch_y_list, index_list_y, batch_kmer_y, spectra_boundaries)

#         # add these these hits to our function scoped set of masses and to the kmer set
#         for k, v in matched_masses_b_batch.items():
#             matched_masses_b[k] += v 
#             kmer_set[k] += batch_kmer_set[k]

#         for k, v in matched_masses_y_batch.items():
#             matched_masses_y[k] += v 
#             kmer_set[k] += batch_kmer_set[k]

#     #update kmers in db to the kmer set
#     db = db._replace(kmers=kmer_set)

#     return (matched_masses_b, matched_masses_y, db)


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

    # turn proteins into namedtuples with .name and .sequence attributes
    basic_prots = [
        BasicProtein(k, v.sequence) for k, v in db.proteins.items()
    ]

    # calculate the max length peptide
    longest_predicted_pep = predicted_len(boundaries[-1][1])

    # get the boundary -> kmer mappings for b and y ions
    matched_masses_b, matched_masses_y, kmer_to_prots = map_spectra_masses.map_boundaries(
        boundaries, 
        basic_prots, 
        longest_predicted_pep
    )

    # replace the db mapping with kmer_to_prots
    db = db._replace(kmers=kmer_to_prots)

    # keep track of the alingment made for every spectrum
    results = {}
    for i, spectrum in enumerate(spectra):
        print(f'\rCreating an alignment for {i+1}/{len(spectra)} [{to_percent(i, len(spectra))}%]', end='')
        results[i] = id_spectrum(spectrum, boundaries, db, mz_mapping, matched_masses_b, matched_masses_y, ppm_tolerance)
        
    return results