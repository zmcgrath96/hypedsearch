from src.objects import Database, Spectrum, Alignments, MPSpectrumID, DEVFallOffEntry
from src.cppModules import gen_spectra
from src.alignment import alignment
from src.utils import ppm_to_da, to_percent, overlap_intervals, hashable_boundaries, is_json, is_file
from src import utils
from src.scoring import scoring, mass_comparisons
from src.preprocessing import digestion, merge_search, preprocessing_utils
from src import database
from src.file_io import JSON

import time
import multiprocessing as mp
import copy
import json

# top results to keep for creating an alignment
TOP_X = 50

def id_spectrum(
    spectrum: Spectrum, 
    db: Database,
    b_hits: dict, 
    y_hits: dict,
    ppm_tolerance: int, 
    precursor_tolerance: int, 
    truth=None, 
    fall_off=None
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
        precursor_tolerance:(int) the toleraence to allow when matching precursor
    Outputs:
        (Alignments) the created alignments for this spectrum
    '''
    # convert the ppm tolerance of the precursor to an int for the rest of the time
    precursor_tolerance = utils.ppm_to_da(spectrum.precursor_mass, precursor_tolerance)

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

    # if fall off and truth are not none, check to see that we can still make the truth seq
    if truth is not None and fall_off is not None:

        # pull out id, hybrid, and truth seq to make it easier
        _id = spectrum.id
        truth_seq = truth[_id]['sequence']
        is_hybrid = truth[_id]['hybrid']

        if not utils.DEV_contains_truth_parts(truth_seq, is_hybrid, filtered_b, filtered_y):

            # add some metadata about what we kept and what fell off
            metadata = {
                'top_x_b_hits': filtered_b, 
                'top_x_y_hits': filtered_y, 
                'excluded_b_hits': [x[0] for x in b_results[keep_b_count:]],
                'excluded_y_hits': [x[0] for x in y_results[keep_y_count:]], 
                'cut_off_b_score': b_results[keep_b_count - 1][1], 
                'cut_off_y_score': y_results[keep_y_count - 1][1]
            }

            # make dev fall off object and add to fall off
            fall_off[_id] = DEVFallOffEntry(
                is_hybrid, 
                truth_seq, 
                'top_x_filtering', 
                metadata
            )

            # skip this entry all together
            return Alignments(spectrum, [])

    # create an alignment for the spectrum
    return alignment.attempt_alignment(
        spectrum, 
        db, 
        filtered_b, 
        filtered_y, 
        ppm_tolerance=ppm_tolerance, 
        precursor_tolerance=precursor_tolerance,
        n=5, 
        truth=truth, 
        fall_off=fall_off
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
    precursor_tolerance=10, 
    digest='',
    missed_cleavages=0,
    cores=1,
    DEBUG=False, 
    truth_set='', 
    output_dir=''
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
        precursor_tolerance:    (int) the tolerance to allow when matching precusor masses. Default=10
        cores:                  (int) the number of cores allowed to use
    Outputs:
        dict containing the results. 
        All information is keyed by the spectrum file name with scan number appended 
        and the values are list of boundariesequenceAligment objects
    '''
    DEV = False
    truth = None

    # for dev use only. If a truth set is passed in, we can check where results
    # drop off. 
    if is_json(truth_set) and is_file(truth_set):
        DEV = True
        print(
            '''
DEV set to True. 
Tracking when correct answer falls off. 
Results are stored in a json named 'fall_off.json' in the specified output directory
File will be of the form

    {
        spectrum_id: {
            hybrid: bool, 
            truth_sequence: str, 
            fall_off_operation: str, 
        }
    }
            '''
        )
        # load in the truth set
        truth = json.load(open(truth_set, 'r'))

    fall_off = None

    # build/load the database
    verbose and print('Loading database...')
    db = database.build(database_file)
    verbose and print('Done')

    
    # load all of the spectra
    verbose and print('Loading spectra...')
    spectra, boundaries, mz_mapping = preprocessing_utils.load_spectra(
        spectra_files, 
        ppm_tolerance,
        peak_filter=peak_filter, 
        relative_abundance_filter=relative_abundance_filter
    )
    verbose and print('Done')

    # get the boundary -> kmer mappings for b and y ions
    matched_masses_b, matched_masses_y, db = merge_search.match_masses(boundaries, db, max_peptide_len)

    # keep track of the alingment made for every spectrum
    results = {}

    if DEV:
        fall_off = {}
        fall_off = mp.Manager().dict()
        truth = mp.Manager().dict(truth)

    # if we only get 1 core, don't do the multiporcessing bit
    if cores == 1:
        # go through and id all spectra
        for i, spectrum in enumerate(spectra):

            print(f'Creating alignment for spectrum {i+1}/{len(spectra)} [{to_percent(i+1, len(spectra))}%]', end='\r')

            # get b and y hits
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

            # pass it into id_spectrum
            results[spectrum.id] = id_spectrum(
                spectrum, 
                db, 
                b_hits, 
                y_hits, 
                ppm_tolerance, 
                precursor_tolerance,
                truth, 
                fall_off
            )

    else:

        print('Initializing other processors...')
        results = mp.Manager().dict()

        if DEV:
            fall_off = mp.Manager().dict()
            truth = mp.Manager().dict(truth)

        # start up processes and queue for parallelizing things
        q = mp.Manager().Queue()
        num_processes = cores
        ps = [
            mp.Process(
                target=mp_id_spectrum, 
                args=(q, copy.deepcopy(db)._replace(b_hits={}, y_hits={}), results, fall_off, truth)
            ) for _ in range(num_processes) 
        ]

        # start each of the process
        for p in ps:
            p.start()
        print('Done.')

        # go through and id all spectra
        for i, spectrum in enumerate(spectra):
            # get b and y hits
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

            # create a named tuple to put in the database
            o = MPSpectrumID(b_hits, y_hits, spectrum, i, ppm_tolerance, precursor_tolerance, 5)
            q.put(o)

        while len(results) < len(spectra):
            print(f'\rCreating an alignment for {len(results)}/{len(spectra)} [{to_percent(len(results), len(spectra))}%]', end='')
            time.sleep(1)

        # now send 'exit' message to all our processes
        [q.put('exit') for _ in range(num_processes)]

        # join them
        for p in ps:
            p.join()

    # if we have set DEV, we need to dump this to a json
    if DEV:
        output_dir = output_dir + '/' if output_dir[-1] != '/' else output_dir

        safe_write_fall_off = {}

        # we need to convert all our DEVFallOffEntries to dicts
        for k, v in fall_off.items():
            safe_write_fall_off[k] = v._asdict()

        JSON.save_dict(output_dir + 'fall_off.json', safe_write_fall_off)
        
    return results

def mp_id_spectrum(
    input_q: mp.Queue, 
    reduced_db_cp: Database, 
    results: dict, 
    fall_off=None, 
    truth=None
    ) -> None:
    '''
    Multiprocessing function for id a spectrum. Each entry in the 
    input_q must be an object of some type with values
        b_hits: list of kmers
        y_hits: list of kmers
        spectrum: Spectrum
        spectrum_id: int
        db: Database
        ppm_tolerance: int
        n: int

    Stores all results in the results dict
    '''
    while True:

        # wait to get something from the input queue
        next_entry = input_q.get(True)

        # if it says 'exit', quit
        if next_entry == 'exit':
            return 

        # if fall off is not none, see if we have the correct value in here
        if truth is not None and fall_off is not None:

            # pull out the id to make it easier
            _id = next_entry.spectrum.id 

            # pull out the truth sequence and the hybrid bool 
            truth_seq = truth[_id]['sequence']
            is_hybrid = truth[_id]['hybrid']

            # see if we still have the correct results
            if not utils.DEV_contains_truth_parts(truth_seq, is_hybrid, next_entry.b_hits, next_entry.y_hits):

                # add some metadata. Add the b and y hits we DID have
                metadata = {
                    'initial_b_candidates': next_entry.b_hits, 
                    'initial_y_candidates': next_entry.y_hits
                }
                
                # create the fall off dev object
                fall_off[_id] = DEVFallOffEntry(
                    is_hybrid, truth_seq, 'mass_matching', metadata
                )
            
                # add an empty entry to the results
                results[_id] = Alignments(next_entry.spectrum, [])
                continue

        # otherwise run id spectrum 
        results[next_entry.spectrum_id] = id_spectrum(
            next_entry.spectrum, 
            reduced_db_cp, 
            next_entry.b_hits, 
            next_entry.y_hits, 
            next_entry.ppm_tolerance, 
            next_entry.precursor_tolerance,
            truth, 
            fall_off
        )