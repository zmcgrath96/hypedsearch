from src.objects import Database, Spectrum, Alignments, MPSpectrumID, DEVFallOffEntry
from src.cppModules import gen_spectra
from src.alignment import alignment
from src.utils import ppm_to_da, to_percent, overlap_intervals, hashable_boundaries, is_json, is_file
from src import utils
from src.scoring import scoring, mass_comparisons
from src.preprocessing import merge_search, preprocessing_utils
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
    n: int,
    digest_type: str = '',
    truth: dict = None, 
    fall_off: dict = None, 
    is_last: bool = False
    ) -> Alignments:
    '''Given the spectrum and initial hits, start the alignment process for 
    the input spectrum

    :param spectrum: observed spectrum in question
    :type spectrum: Spectrum
    :param db: Holds all the source sequences
    :type db: Database
    :param b_hits: all k-mers found from the b-ion search
    :type b_hits: list
    :param y_hits: all k-mers found from the y-ion search
    :type y_hits: list
    :param ppm_tolerance: the parts per million error allowed when trying to match masses
    :type ppm_tolerance: int
    :param precursor_tolerance: the parts per million error allowed when trying to match
        precursor masses
    :type percursor_tolerance: int
    :param n: the number of alignments to save
    :type n: int
    :param digest_type: the digest performed on the sample
        (default is '')
    :type digest_type: str
    :param truth: a set of id keyed spectra with the desired spectra. A better description of what this looks like can be 
        seen in the param.py file. If left None, the program will continue normally
        (default is None)
    :type truth: dict
    :param fall_off: only works if the truth param is set to a dictionary. This is a dictionary (if using multiprocessing, 
        needs to be process safe) where, if a sequence loses the desired sequence, a key value pair of spectrum id, 
        DevFallOffEntry object are added to it. 
        (default is None)
    :type fall_off: dict
    :param is_last: Only works if DEV is set to true in params. If set to true, timing evaluations are done. 
        (default is False)
    :type is_last: bool

    :returns: Alignments for the spectrum. If no alignment can be created, and empty Alignments object is inserted
    :rtype: Alignments
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

    # count the number of kmers that have the highest value
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
        n=n, 
        truth=truth, 
        fall_off=fall_off, 
        is_last=is_last
    )


def id_spectra(
    spectra_files: list, 
    database_file: str, 
    verbose: bool = True, 
    min_peptide_len: int = 5, 
    max_peptide_len: int = 20, 
    peak_filter: int = 0, 
    relative_abundance_filter: float = 0.0,
    ppm_tolerance: int = 20, 
    precursor_tolerance: int = 10, 
    digest: str = '',
    cores: int = 1,
    n: int = 5,
    DEBUG: bool = False, 
    truth_set: str = '', 
    output_dir: str = ''
) -> dict:
    '''Load in all the spectra and try to create an alignment for every spectrum

    :param spectra_files: file names of input spectra
    :type spectra_files: list
    :param database_file: file name of the fasta database
    :type database_file: str
    :param verbose: print progress to the console. 
        (default is True)
    :type verbose: bool
    :param min_peptide_len: the minimum length alignment to create
        (default is 5)
    :type min_peptide_len: int
    :param max_peptide_len: the maximum length alignment to create
        (default is 20)
    :type max_peptide_len: int
    :param peak_filter: If set to a number, this metric is used over the relative abundance filter. 
        The most abundanct X peaks to use in the alignment. 
        (default is 0)
    :type peak_filter: int
    :param relative_abundance_filter: If peak_filter is set, this parameter is ignored. The 
        relative abundance threshold (in percent as a decimal) a peak must be of the total 
        intensity to be used in the alignment. 
        (default is 0.0)
    :type relative_abundance_filter: float
    :param ppm_tolerance: the parts per million error allowed when trying to match masses
        (default is 20)
    :type ppm_tolerance: int
    :param precursor_tolerance: the parts per million error allowed when trying to match
        a calculated precursor mass to the observed precursor mass
        (default is 10)
    :type precurosor_tolerance: int
    :param digest: the type of digest used in the sample preparation. If left blank, 
        a digest-free search is performed. 
        (default is '')
    :type digest: str
    :param cores: the number of cores allowed to use in running the program. If a number 
        provided is greater than the number of cores available, the maximum number of 
        cores is used. 
        (default is 1)
    :type cores: int
    :param n: the number of aligments to keep per spectrum. 
        (default is 5)
    :type n: int
    :param DEBUG: DEVELOPMENT USE ONLY. Used only for timing of modules. 
        (default is False)
    :type DEBUG: bool
    :param truth_set: the path to a json file of the desired alignments to make for each spectrum. 
        The format of the file is {spectrum_id: {'sequence': str, 'hybrid': bool, 'parent': str}}. 
        If left an empty string, the program proceeds as normal. Otherwise results of the analysis
        will be saved in the file 'fall_off.json' saved in the output directory specified.
        (default is '')
    :type truth_set: str
    :param output_dir: the full path to the output directory to save all output files.
        (default is '')
    :type output_dir: str

    :returns: alignments for all spectra save in the form {spectrum.id: Alignments}
    :rtype: dict
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

    # if we only get 1 core, don't do the multiprocessing bit
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

            is_last = DEBUG and i == len(spectra) - 1

            # pass it into id_spectrum
            results[spectrum.id] = id_spectrum(
                spectrum, 
                db, 
                b_hits, 
                y_hits, 
                ppm_tolerance, 
                precursor_tolerance,
                n,
                digest_type=digest,
                truth=truth, 
                fall_off=fall_off, 
                is_last=is_last
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
                args=(q, copy.deepcopy(db), results, fall_off, truth)
            ) for _ in range(num_processes) 
        ]

        # start each of the process
        for p in ps:
            p.start()
        print('Done.')

        # go through and id all spectra
        for i, spectrum in enumerate(spectra):
            print(f'\rStarting job for {i+1}/{len(spectra)} [{to_percent(i+1, len(spectra))}%]', end='')
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
            o = MPSpectrumID(
                b_hits, 
                y_hits, 
                spectrum, 
                ppm_tolerance, 
                precursor_tolerance, 
                n, 
                digest
            )
            
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
    db_copy: Database, 
    results: dict, 
    fall_off: dict = None, 
    truth: dict = None
    ) -> None:
    '''Multiprocessing function for to identify a spectrum. Each entry in the 
    input_q must be a MPSpectrumID object

    :param input_q: a queue to pull MPSpectrumID objects from for analysis
    :type input_q: mp.Queue
    :param db_copy: a copy of the original database for alignments
    :type db_copy: Database
    :param results: a multiprocesses safe dictionary to save the alignments in
    :type results: dict
    :param truth_set: dictionary containing all the desired alignments to make. 
        The format of the file is {spectrum_id: {'sequence': str, 'hybrid': bool, 'parent': str}}. 
        If left as None, the program will continue as normal
        (default is None)
    :type truth_set: dict
    :param fall_off: only used if the truth_set param is set to a valid json. Must be a multiprocess
        safe dictionary to store the fall off information to
    :type fall_off: dict

    :returns: None
    :rtype: None
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
        results[next_entry.spectrum.id] = id_spectrum(
            next_entry.spectrum, 
            db_copy, 
            next_entry.b_hits, 
            next_entry.y_hits, 
            next_entry.ppm_tolerance, 
            next_entry.precursor_tolerance,
            next_entry.n,
            truth, 
            fall_off
        )