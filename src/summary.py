import json
from src.utils import make_dir, make_valid_dir_string
from src.file_io import JSON

SUMMARY_NAME = 'summary'

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