from src.file_io import mzML
from src.objects import Database, Spectrum, Alignments
from src.sequence.gen_spectra import gen_spectrum, gen_min_ordering
from src.alignment import alignment
from src.utils import ppm_to_da, to_percent
from src.scoring import scoring, mass_comparisons
from src import database

from math import ceil
from collections import defaultdict

import bisect

# top results to keep for creating an alignment
TOP_X = 10

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
    b_hits = defaultdict(set)
    y_hits = defaultdict(set)

    # the maximum lenth peptide we could possibly find. Calculation is:
    # the ceiling of ( the maximum observed mass / the smallest singly charged mass (G))
    max_len = int(ceil(reduced_spectra[-1]/ 57.021464))

    def find_kmer_hits(kmer: str, prot_name: str) -> None:
        '''
        Internal function for adding both ions and charges to the tree 
        and to the dictionaries too keep track of sequences
        '''
        for ion in 'by':
            for charge in [1, 2]:
   
                spec = gen_spectrum(kmer, ion=ion, charge=charge)['spectrum']

                for c, mass in enumerate(spec):

                    # calculate the upper and lower bounds in the search
                    da_tol = ppm_to_da(mass, ppm_tol)
                    lb = mass - da_tol
                    ub = mass + da_tol

                    # get the index of the value to the left of any existing values of 
                    # the lower bound. Any value after this postion is greater or equal to lower bound
                    beginning_entry = bisect.bisect_left(reduced_spectra, lb)

                    if kmer[-c-1:] == 'EAPNFEANTTIGRIRFH' and ion == 'y':
                        print('potato')

                    # see if the NEXT value is in the range. If so, keep the kmer
                    to_insert = False
                    if beginning_entry + 1 < len(reduced_spectra) and reduced_spectra[beginning_entry] <= ub:
                        to_insert = True
                        offset = 0
                        while reduced_spectra[beginning_entry + offset] <= ub:
                            
                            # increment the hits and the longest kmer found
                            if ion == 'b':
                                b_hits[reduced_spectra[beginning_entry + offset]].add(kmer[:c+1])
                            else:
                                y_hits[reduced_spectra[beginning_entry + offset]].add(kmer[-c-1:])

                            offset += 1

                    if to_insert:
                        if ion == 'b':
                            db.tree.insert(prot_name, kmer[:c+1])
                        else:
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
    db = db._replace(b_hits=b_hits)
    db = db._replace(y_hits=y_hits)

    print('Done.')
    return db


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
            x for spectrum in these_spectra for x in spectrum.spectrum
        ]))

    return (all_spectra, sorted(list(set(linear_spectra))))


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
    db = reduce_search_space(db, reduced_spectra, ppm_tolerance)

    # keep track of each of the results
    b_results = defaultdict(list)
    y_results = defaultdict(list)

    # go through each spectrum, sort their results, and take the top X hits to try and align
    results = {}
    for i, spectrum in enumerate(spectra):
        print(f'\rCreating an alignment for {i+1}/{len(spectra)} [{to_percent(i, len(spectra))}%]', end='')
        # b_results[i].sort(key=lambda x: x[1], reverse=True)
        # y_results[i].sort(key=lambda x: x[1], reverse=True)
        b_hits, y_hits = [], []
        for mz in spectrum.spectrum:
            if mz in db.b_hits:
                b_hits += list(db.b_hits[mz])
            if mz in db.y_hits:
                y_hits += list(db.y_hits[mz])

        # remove any duplicates
        b_hits = list(set(b_hits))
        y_hits = list(set(y_hits))

        # score and sort these results
        b_results = sorted([
            (
                kmer, 
                mass_comparisons.optimized_compare_masses(spectrum.spectrum, gen_spectrum(kmer, ion='b')['spectrum'])
            ) for kmer in b_hits], 
            key=lambda x: x[1], 
            reverse=True
        )
        y_results = sorted([
            (
                kmer, 
                mass_comparisons.optimized_compare_masses(spectrum.spectrum, gen_spectrum(kmer, ion='y')['spectrum'])
            ) for kmer in y_hits], 
            key=lambda x: x[1], 
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
        results[i] = alignment.attempt_alignment(
            spectrum, 
            db, 
            filtered_b, 
            filtered_y, 
            ppm_tolerance=ppm_tolerance, 
            n=3, 
            scoring_alg='ion'
        )
    return results