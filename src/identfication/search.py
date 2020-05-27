from src.identfication import filtering, kmer_extension
from src.scoring import scoring
from src.types.database import Database
from src.types.objects import Spectrum, Kmer, ScoredKmer, MassSequence
from src.spectra.gen_spectra import gen_spectrum

################### Constants ###################
BASE_K = 3
SDEVS = 4
STALL_LENGTH = 3
#################################################
    

def get_interesting_protein_positions(database: Database, mers: list) -> list:
    '''
    Get the interesting protein entries and positions from a list of mers

    Inputs:
        database:   instance of class Database
        mers:       list of strings kmers found
    Outputs:
        list of Kmer namedtuple instances
    '''
    interesting_spots = []
    for mer in mers:
        mdom = database.get_metadata_of_mer(mer)
        [interesting_spots.append(x) for x in mdom]
    return interesting_spots

def extend_spots_of_interest(spectrum: Spectrum, database: Database, spots_of_interest: list, ion_type: str, n=3) -> list:
    '''
    Look through the spots of interest and extend them/score them as much as possible

    Inputs:
        spectrum:           Spectrum namedtuple instance
        database:           instance of the Database class
        spots_of_interest:  list of tuples of Kmer namedtuple instances
        ion_type:           str the ion type we are scoring. Should be either 'b' or 'y'
    kwargs:
        n:                  int the top number of results to return. Default=3 
    Outputs:
        list of sorted entries of ScoredKmer instances
    '''
    if ion_type.lower() not in ['b', 'y']:
        return {}
    sort_by_key = 'b_score' if ion_type.lower() == 'b' else 'y_score'
    extended = []
    local_entries = {}
    for soi in spots_of_interest:
        # get the entry from the database
        if soi.protein not in local_entries:
            local_entries[soi.protein] = database.get_entry_by_name(soi.protein)
        # extend the soi from either starting or ending depending on the ion
        b_score, y_score = scoring.score_subsequence(spectrum.spectrum, soi.sequence)
        sk = ScoredKmer(b_score, y_score, soi)
        extended_kmer = kmer_extension.extend_kmer(spectrum, local_entries[soi.protein], sk, ion_type, STALL_LENGTH)
        extended.append(extended_kmer)
    extended.sort(key=lambda x: getattr(x, sort_by_key), reverse=True)
    return extended[:n]


def search_database(spectrum: Spectrum, database: Database, base_k_mers: list,  n=3, ppm_tolerance=20) -> (dict, dict):
    '''
    Search through the proteins to find the top n sequences that describe the spectrum

    Inputs:
        spectrum:       (Spectrum) namedtuple instance
        database:       (Database) class instance
        base_k_mers:    (list) list of the top scoring kmers to extend
    kwargs:
        n:              (int) top number of results to return for both b and y ion scores. Default=3
        ppm_tolerance:  (int) mass tolerance to look out for. Default=20
    Outptus:
        (b_dict, y_dict)
        Both of these ion dictionaries have the top n entries keyed by their rank (0 to n-1, 0 being the best)
        with values being ScoredKmer namedtuple instances
    '''   
    # we want to go through ONLY the proteins that have these subsequences in them
    b_spots_of_interest = get_interesting_protein_positions(database, base_k_mers)
    y_spots_of_interest = get_interesting_protein_positions(database, base_k_mers)
    # we have the interesting proteins and associated starting/ending positions
    # try and extend these as much as possible
    best_b = extend_spots_of_interest(spectrum, database, b_spots_of_interest, 'b', n)
    best_y = extend_spots_of_interest(spectrum, database, y_spots_of_interest, 'y', n)
    
    # dont go out of range
    iter_range = min(n, len(best_b), len(best_y))

    indexed_b = {i: best_b[i] for i in range(iter_range)}
    indexed_y = {i: best_y[i] for i in range(iter_range)}
    return indexed_b, indexed_y