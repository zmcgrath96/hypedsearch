from src.alignment import filtering, kmer_extension
from src.scoring import scoring
from src.database.database import Database
from src.modules.spectrum import Spectrum
from src.modules.kmer import Kmer
from src.modules.scores import Scores
from typing import Iterable

################### Constants ###################
BASE_K = 3
SDEVS = 4
STALL_LENGTH = 3
#################################################

def find_interesting_kmers(pepspec: Spectrum, mers: list) -> (list, list):
    '''
    Go through a list of mers and find the ones that are considered interesting

    Inputs:
        pepspec:    Spectrum object
        mers:       list of Kmer objects
    Outputs:
        (b_anchors, y_anchors): list of dicts of the form {b_score: float, y_score: float, mer: str}
    '''
    # for ever mer, score the subsequence
    mer_scores = []
    for mer in mers:
        # create a scores object and score it, appending it to the mer scores
        s = Scores(pepspec, mer)
        s.score()
        mer_scores.append(s)
    b_anchors = filtering.stddev_filter(mer_scores, 'b_score', SDEVS)
    y_anchors = filtering.stddev_filter(mer_scores, 'y_score', SDEVS)
    return b_anchors, y_anchors

def get_interesting_protein_positions(database: Database, mers: Iterable) -> list:
    '''
    Get the interesting protein entries and positions from a list of mers

    Inputs:
        database:   instance of class Database
        mers:       iterable that returns kmers
    Outputs:
        list of tuples with entries (protein name: str, starting_position: int, ending_position: int)
    '''
    interesting_spots = []
    for mer in mers:
        mdom = database.get_metadata_of_mer(mer)
        [interesting_spots.append(x) for x in mdom]
    return interesting_spots

def get_best_spots_of_interest(spectrum: Spectrum, database: Database, spots_of_interest: list, ion_type: str, n=3) -> list:
    '''
    Look through the spots of interest and extend them/score them as much as possible

    Inputs:
        spectrum:           Spectrum object
        database:           instance of the Database class
        spots_of_interest:  list of tuples of the form (protein_name: str, starting_position, ending_position)
        ion_type:           str the ion type we are scoring. Should be either 'b' or 'y'
    kwargs:
        n:                  int the top number of results to return. Default=3 
    Outputs:
        list of Scores objects
    '''
    if ion_type.lower() not in ['b', 'y']:
        return {}
    sort_by_key = 'b_score' if ion_type.lower() == 'b' else 'y_score'
    extended = []
    local_entries = {}
    for soi in spots_of_interest:
        # get the entry from the database
        if soi[0] not in local_entries:
            local_entries[soi[0]] = database.get_entry_by_name(soi[0])
        entry = local_entries[soi[0]]
        kmer = Kmer(entry.sequence[soi[1]:soi[2]], entry.name, soi[1])
        # score it 
        s = Scores(spectrum, kmer)
        s.score()
        # extend the soi from either starting or ending depending on the ion
    
        extended_kmer = kmer_extension.extend_kmer(spectrum, entry.sequence, s, ion_type, STALL_LENGTH)
        extended.append(extended_kmer)
    extended.sort(key=lambda x: getattr(x, sort_by_key), reverse=True)
    return extended[:n]


def search_database(spectrum: Spectrum, database: Database, n=3) -> (dict, dict):
    '''
    Search through the proteins to find the top n sequences that describe the spectrum

    Inputs:
        spectrum:   Spectrum object
        database:   Database class instance
    kwargs:
        n:      int top number of results to return for both b and y ion scores. Default=3
    Outptus:
        (b_dict, y_dict)
        Both of these ion dictionaries have the top n entries keyed by their rank (0 to n-1, 0 being the best)
         {
             0: {starting_position: int, ending_position: int, k: int, b_score: float, y_score: float, sequence: str, protein_name: str}
             ...
             n-1: {...}
         }
    '''
    # go through the metadata in the database to find the interesting kmers
    metadata = database.metadata
    # the metadata keys are the kmers
    mers = [Kmer(mer, entry[0], entry[1]) for mer, entry in metadata]
    # go through all of the mers and find only the interesting ones
    b_anchors, y_anchors = find_interesting_kmers(spectrum, mers)
    # sort them so that we only take the most interesting ones
    b_anchors.sort(key=lambda x: x.b_score, reverse=True)
    y_anchors.sort(key=lambda x: x.y_score, reverse=True)
    # take the kmers from them
    b_anchor_mers = [x.kmer.sequence for x in b_anchors]
    y_anchor_mers = [x.kmer.sequence for x in y_anchors]
    # we have the interesting b anchors and y anchors
    # we want to go through ONLY the proteins that have these subsequences in them
    b_spots_of_interest = get_interesting_protein_positions(database, b_anchor_mers)
    y_spots_of_interest = get_interesting_protein_positions(database, y_anchor_mers)
    # we have the interesting proteins and associated starting/ending positions
    # try and extend these as much as possible
    best_b = get_best_spots_of_interest(spectrum ,database, b_spots_of_interest, 'b', n)
    best_y = get_best_spots_of_interest(spectrum, database, y_spots_of_interest, 'y', n)
    
    # dont go out of range
    iter_range = min(n, len(best_b), len(best_y))

    indexed_b = {i: best_b[i] for i in range(iter_range)}
    indexed_y = {i: best_y[i] for i in range(iter_range)}
    return indexed_b, indexed_y