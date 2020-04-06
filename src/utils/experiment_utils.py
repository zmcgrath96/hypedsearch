import copy

def split_exp_by_ion(exp: dict, ion: str) -> dict:
    '''
    Make a copy of the experiment split by ion information

    Inputs:
        exp:    dictionary with expeiment summary information
        ion:    string ion type. Possible types are {'b', 'y'}
    Outpus:
        dict same structure as experiment just without other ion information
    '''
    ionized_exp = {}
    ionized_exp['experiment_info'] = copy.deepcopy(exp['experiment_info'])
    ionized_exp['experiment'] = {}
    for pep_name, pep in exp['experiment'].items():
        ionized_exp['experiment'][pep_name] = {}
        for prot_name, prot in pep.items():
            if prot_name == 'analysis':
                ionized_exp['experiment'][pep_name][prot_name] = copy.deepcopy(prot)
                ionized_exp['experiment'][pep_name][prot_name]['ranks']['ranks'] = ionized_exp['experiment'][pep_name][prot_name]['ranks']['ranks'][ion]
                ionized_exp['experiment'][pep_name][prot_name]['sequence_predictions'] = ionized_exp['experiment'][pep_name][prot_name]['sequence_predictions'][ion]
                continue
            ionized_exp['experiment'][pep_name][prot_name] = {}
            for k, scores in prot.items():
                ionized_exp['experiment'][pep_name][prot_name][k] = scores[ion] 

    return ionized_exp

def experiment_has_ion_types(exp: dict) -> bool:
    '''
    Determines if the experiment split scores by ion type

    Inputs:
        exp:    dictionary containing experiment information
    Outputs:
        True if ion types detected False otherwise
    '''
    first_pep = exp['experiment'][next(iter(exp['experiment']))]
    return set(['b', 'y']) & set(first_pep['analysis']['ranks']['ranks'].keys())