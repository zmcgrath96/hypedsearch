from src.utils import make_dir, make_valid_dir_string
from src.file_io import JSON
from src.types.objects import Alignments

import pandas as pd
import json

SUMMARY_NAME = 'summary'
HYBRID_PREFIX = 'hybrid_'

def json_file(results: dict, output_dir='./') -> None:
    '''
    Generate a summary json file for the results made

    Inputs:
        results: (dict) containing the results made. The key of each entry is the name
                        and each entry should be an Alignments namedtuple
    kwargs:
        output_dir: (str) path to the output directory. Default=./
    Outputs:
        None
    '''
    json_file_name = output_dir + SUMMARY_NAME 
    dictified = {}
    for name, alignment in results.items():
        dictified[name] = {
            'spectrum': alignment.spectrum._asdict(), 
            'alignments': [x._asdict() for x in alignment.alignments]
        }
    JSON.save_dict(json_file_name, dictified)

def tsv_file(results: dict, output_dir='./') -> None:
    '''
    Write the results of the experiment to 2 tsv files. One tsv file is for 
    non hybrid identified sequences, the other is for hybrid identified sequences.

    Inputs: 
        results:    (dict) results of the search. The key of each entry is the name
                            of each entry and the value is an Alignments namedtuple 
    kwargs:
        output_dir: (str) path to the directory to save the tsvs. Default = ./
    Outputs:
        None
    '''
    mac = 0

    # seperate the hybrids from the nonhybrids
    hybrids, nonhybrids = [], []
    alignment: Alignments
    for name, alignment in results.items():
        if len(alignment.alignments) == 0:
            mac += 1
            continue

        topalignment = alignment.alignments[0]._asdict()
        topalignment['entry name'] = name
        topalignment['id'] = alignment.spectrum.id

        if 'hybrid_sequence' in topalignment:
            hybrids.append(topalignment)
        else:
            nonhybrids.append(topalignment)

    # move to pandas dataframe for easy writing
    hybridresults = pd.DataFrame(hybrids)
    with open(f'{output_dir + HYBRID_PREFIX + SUMMARY_NAME}.tsv', 'w') as ho:
        ho.write(hybridresults.to_csv(sep='\t'))

    del hybridresults
    del hybrids

    nonhybridresults = pd.DataFrame(nonhybrids)
    with open(f'{output_dir + SUMMARY_NAME}.tsv', 'w') as nho:
        nho.write(nonhybridresults.to_csv(sep='\t'))

    print(f'Could not make an alignment for {mac}/{len(results)} spectra ({int(100 * mac / len(results))}%)')

def generate(alignments: dict, output_dir='./') -> None:
    '''
    Generate a summary text and json file for the alignments made

    Inputs:
        alignments: dict containing the alignments made. The key of each entry is the name
                    of the file appended with scan number, and the values should be Alignments
    kwargs:
        output_dir: str path to the output directory. Default=./
    Outputs:
        None
    '''
    output_dir = make_valid_dir_string(output_dir)
    make_dir(output_dir)

    json_file(alignments, output_dir)
    tsv_file(alignments, output_dir=output_dir)