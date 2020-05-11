import json
from src.utils.utils import make_dir, make_valid_dir_string
from src.file_io import JSON

SUMMARY_NAME = 'summary'

def json_file(alignments: dict, output_dir='./') -> None:
    '''
    Generate a summary json file for the alignments made

    Inputs:
        alignments: dict containing the alignments made. The key of each entry is the name
                    and each entry should contain the following:
                        {spectrum: list, alignments: list, scan_no: int}
    kwargs:
        output_dir: str path to the output directory. Default=./
    Outputs:
        None
    '''
    json_file_name = output_dir + SUMMARY_NAME 
    JSON.save_dict(json_file_name, alignments)

def text_file(alignments: dict, output_dir='./') -> None:
    '''
    Generate a summary text file for the alignments made

    Inputs:
        alignments: dict containing the alignments made. The key of each entry is the name
                    and each entry should contain the following:
                        {spectrum: list, alignments: list, scan_no: int}
    kwargs:
        output_dir: str path to the output directory. Default=./
    Outputs:
        None
    '''
    txt_file_name = output_dir + SUMMARY_NAME + '.txt'
    entries = ['filename', 'scan_no', 'protein_name', 'sequence', 'length', 'starting_position', 'ending_position', 'b_score', 'y_score', 'confidence']
    line_template = ','.join(['{}' for _ in range(10)])
    header_row = line_template.format(*entries)

    o = open(txt_file_name, 'w')
    o.write(header_row)

    for name, spectrum in alignments.items():
        # each alignmetn has the following structure
        # {starting_position: int, ending_position: int, sequence: str, length: int, b_score: float, y_score: float, confidence: float, protein_name: string, spectrum: list}
        filename = name.split('/')[-1].split('_')[0]
        als = spectrum['alignments']
        als.sort(key=lambda x: x['confidence'], reverse=True)
        best_alignment = als[0]
        entry = line_template.format(filename, spectrum['scan_no'], *[best_alignment[x] for x in entries[2:]])
        o.write(entry)

def generate(alignments: dict, output_dir='./') -> None:
    '''
    Generate a summary text and json file for the alignments made

    Inputs:
        alignments: dict containing the alignments made. The key of each entry is the name
                    and each entry should contain the following:
                        {spectrum: list, alignments: list, scan_no: int}
    kwargs:
        output_dir: str path to the output directory. Default=./
    Outputs:
        None
    '''
    output_dir = make_valid_dir_string(output_dir)
    make_dir(output_dir)

    # text_file(alignments, output_dir)
    json_file(alignments, output_dir)