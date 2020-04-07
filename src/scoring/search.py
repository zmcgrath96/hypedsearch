from file_io import fasta, csv
from scoring.comparisons import compare_sequence_spectra, compare_spectra_sequence_ion_type

def ion_matches(spectrum: list, database: list, n=3, ion='b') -> dict:
    '''
    Find the top n matches for each ion type

    Inputs:
        spectrum:   list of floats the mass spectrum
        database:   list of dictionaries from a fasta file of the form
                    {
                        'name': str, 'sequence': str, 'identifier': str
                    }
    Outputs:
        dictionary of the form:
            {
                0: {protein_name: str, peptide_sequence: str, score: float, starting_position: int, ending_position: int}, 
                ...
                n-1: {...}
            }
    '''
    top_n_matches = []


def search_proteins(spectrum: dict, database: list) -> (dict, dict):
    '''
    MARCH 11 202
    Find the highest scoring protein from a spectrum in a database
    Uses both b and y ion scores to find it. 

    Inputs:
        spectrum:   dictionary of the form {'level': , 'scan_no': , 'spectrum': }
        database:   list of dictionaries read from a fasta file of form
                    {
                        'name': str, 'sequence': str, 'identifier': str
                    }
    Outputs:
        (b_match, y_match)
        both of these matches have the structure: 
        {
            'protein_id': any,
            'protein_name': str,
            'protein_sequence': str,
            'ms_level': int,
            'scan_no': int,
            'score': float
        }
    '''
    best_match_b = {}
    best_match_y = {}
    for prot in database:
        b_score = compare_spectra_sequence_ion_type(spectrum['spectrum'], prot['sequence'], 'b')
        y_score = compare_spectra_sequence_ion_type(spectrum['spectrum'], prot['sequence'], 'y')
        if best_match_b == {} or best_match_b['score'] < b_score:
            best_match_b = {
                'protein_id': prot['identifier'],
                'protein_name': prot['name'], 
                'protein_sequence': prot['sequence'], 
                'ms_level': spectrum['level'],
                'scan_no': spectrum['scan_no'], 
                'score': b_score
            }
        if best_match_y == {} or best_match_y['score'] < y_score:
            best_match_y = {
                'protein_id': prot['identifier'],
                'protein_name': prot['name'], 
                'protein_sequence': prot['sequence'], 
                'ms_level': spectrum['level'],
                'scan_no': spectrum['scan_no'], 
                'score': y_score
            }

    best_matches = {
        'protein_id_b': best_match_b['protein_id'],
        'protein_name_b': best_match_b['protein_name'],
        'ms_level_b': best_match_b['ms_level'],
        'scan_no_b': best_match_b['scan_no'],
        'score_b': best_match_b['score'],
        'protein_id_y': best_match_y['protein_id'],
        'protein_name_y': best_match_y['protein_name'],
        'ms_level_y': best_match_y['ms_level'],
        'scan_no_y': best_match_y['scan_no'],
        'score_y': best_match_y['score']
    }
    return best_matches
